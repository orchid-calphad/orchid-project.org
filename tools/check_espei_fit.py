#!/usr/bin/env python3
"""Read-only preflight of an explicit fitting selection; requires pycalphad and ESPEI.
Usage: python check_fit.py --bundle espei --selection selection.json
All paths in selection.json are relative to that file, except dataset path (bundle-relative).
No files, weights, phase names or observations are rewritten.
"""
import argparse, collections, hashlib, json, math
from pathlib import Path

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x): return json.dumps(x, sort_keys=True, separators=(',', ':'))
def positive(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and x > 0

def evidence(d):
    o=d.get('orchid', {}); direct=set(); derived=set()
    for ob in o.get('observations', []):
        direct.add(ob['id'])
        for alias in ob.get('also_reported_as', [])+ob.get('also_read_as', []):
            direct.add(alias['id'] if isinstance(alias, dict) else alias)
    for m in o.get('mapping', []):
        if m.get('id'):direct.add(str(m['id']) if str(m['id']).startswith('ORC-') else str(o.get('dft_source'))+':'+str(m['id']))
        direct.update(m.get('also_read_as', []))
    direct.update(o.get('record_ids', []))
    derived.update(p['source_observation'] for p in o.get('points', []) if p.get('source_observation'))
    return direct, derived

def check(bundle, selection):
    from pycalphad import Database
    from espei.datasets import load_datasets
    bundle=bundle.resolve();selection=selection.resolve();manifest=json.loads(selection.read_text());base=selection.parent
    errors=[];reviews=[];notes=[];inputs=[];seen_paths=set();entries=[]
    def err(code, **kw):errors.append(dict(code=code, **kw))
    target=manifest['target'];tp=(base/target['path']).resolve();th=digest(tp)
    if target.get('sha256')!=th:err('target_hash_mismatch')
    db=Database(str(tp));registry={};identities={}
    for p in (bundle.parent/'phases/systems').glob('*.json'):
        for phase in json.loads(p.read_text())['phases']:
            for m in phase.get('models', []):
                registry[(m['tdb'],m['phase_name'])]=phase['id']
                if m.get('sha256')==th:identities[m['phase_name']]=phase['id']
    # Explicit identity assertions for new models remain visible in the receipt.
    assertions=target.get('phase_identity', {})
    for phase, assertion in assertions.items():
        if not isinstance(assertion,dict) or not assertion.get('registry_id') or not str(assertion.get('reason','')).strip():
            err('invalid_phase_identity_assertion',phase=phase);continue
        if phase not in db.phases:err('unknown_asserted_phase',phase=phase);continue
        if phase in identities and identities[phase]!=assertion['registry_id']:err('contradictory_phase_identity',phase=phase)
        else:identities[phase]=assertion['registry_id']
    index={r['path']:r for r in json.loads((bundle/'index.json').read_text())['datasets']}
    if not manifest.get('datasets'):err('empty_selection')
    for item in manifest.get('datasets', []):
        path=item['path'];p=(bundle/path).resolve()
        if not p.is_relative_to(bundle):err('dataset_outside_bundle',path=path);continue
        if path in seen_paths:err('repeated_dataset_path',path=path);continue
        seen_paths.add(path)
        if path not in index:err('dataset_not_indexed',path=path);continue
        raw=json.loads(p.read_text());sha=digest(p)
        if item.get('sha256')!=sha:err('dataset_hash_mismatch',path=path)
        fp=(base/item['fit_file']).resolve() if item.get('fit_file') else p
        d=json.loads(fp.read_text());inputs.append({'path':path,'bundle_sha256':sha,'fit_file':str(fp),'fit_sha256':digest(fp)})
        if {k:v for k,v in d.items() if k!='weight'}!={k:v for k,v in raw.items() if k!='weight'}:
            err('fit_file_changed_beyond_weight',path=path)
        w=d.get('weight',1);weights=w if isinstance(w,list) else [w]
        if not weights or not all(positive(x) for x in weights):err('invalid_weight',path=path)
        override=item.get('weight_override')
        if w!=raw.get('weight',1) or override is not None:
            if not isinstance(override,dict) or override.get('value')!=w or not str(override.get('reason','')).strip():
                err('undocumented_weight_override',path=path)
        output=d['output'];u=raw.get('orchid_uncertainty',{})
        if output!='ZPF':
            uv=u.get('values',[]);bw=raw.get('weight',1);bw=bw if isinstance(bw,list) else [bw]*len(uv)
            if not uv or len(bw)!=len(uv):err('uncertainty_alignment',path=path)
            else:
                for v,weight in zip(uv,bw):
                    sigma=v.get('sigma_effective_used',v.get('budget',{}).get('sigma_effective',v.get('sigma')));default=u.get('default_sigma',{}).get('value')
                    if not positive(sigma) or not positive(default) or not positive(weight) or abs(weight-round(default/sigma,6))>0.00000051:
                        err('invalid_exported_uncertainty_weight',path=path);break
        if not set(d['components'])<=set(db.elements):err('target_components_missing',path=path)
        phases=set(d['phases']);phases.update(d.get('reference_state',{}).get('phases',[]))
        refs=d.get('reference_states',[])
        refs=refs.values() if isinstance(refs,dict) else refs
        for r in refs:
            if isinstance(r,dict) and r.get('phase'):phases.add(r['phase'])
        missing=phases-set(db.phases)
        if missing:err('target_phases_missing',path=path,phases=sorted(missing))
        source=index[path].get('tdb') or Path(path).name.split('.')[0]+'.tdb'
        for phase in d['phases']:
            a=registry.get((source,phase));b=identities.get(phase)
            if a and b and a!=b:err('registry_phase_mismatch',path=path,phase=phase)
            elif not a or not b:reviews.append({'code':'unverified_phase_identity','paths':[path],'phase':phase})
        solver=d.get('solver')
        if solver and not missing:
            if len(d['phases'])!=1:err('solver_requires_one_phase',path=path)
            else:
                ph=db.phases[d['phases'][0]];ratios=solver.get('sublattice_site_ratios',[])
                if len(ratios)!=len(ph.sublattices) or any(not math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-7) for a,b in zip(ratios,ph.sublattices)):
                    err('solver_sublattice_mismatch',path=path)
                for config in solver.get('sublattice_configurations',[]):
                    if len(config)!=len(ph.constituents):err('solver_configuration_mismatch',path=path);continue
                    for species,allowed in zip(config,ph.constituents):
                        species=species if isinstance(species,list) else [species]
                        if not set(species)<={s.name for s in allowed}:err('solver_species_mismatch',path=path)
        try:load_datasets([str(fp)])
        except Exception as exc:err('espei_load_error',path=path,message=str(exc))
        direct,derived=evidence(d)
        payload=canonical({k:d.get(k) for k in ['components','phases','conditions','output','values','reference_state','reference_states','solver']})
        entries.append((path,d,direct,derived,payload))
    for i,(path,d,ids,derived,payload) in enumerate(entries):
        for other,e,other_ids,other_derived,other_payload in entries[:i]:
            pair=sorted([path,other]);overlap=ids&other_ids
            if overlap:err('repeated_evidence',paths=pair,ids=sorted(overlap));continue
            if (derived&(other_ids|other_derived)) or (other_derived&ids):
                reviews.append({'code':'shared_derived_evidence','paths':pair})
            if payload==other_payload:reviews.append({'code':'identical_numeric_constraints','paths':pair})
            od=d.get('orchid',{});oe=e.get('orchid',{})
            if od.get('diagram') and od.get('diagram')==oe.get('diagram'):
                reviews.append({'code':'shared_assessed_diagram','paths':pair})
    acknowledged=[];pending=[]
    for r in reviews:
        matches=[a for a in manifest.get('acknowledgements',[]) if a.get('code')==r['code'] and sorted(a.get('paths',[]))==sorted(r['paths']) and str(a.get('reason','')).strip()]
        (acknowledged if matches else pending).append(dict(r,acknowledgement=matches[0] if matches else None))
    provenance=None;prov=bundle.parent/'phases/tdb_provenance.json'
    if prov.exists():provenance=next((r for r in json.loads(prov.read_text())['records'] if r['sha256']==th),None)
    return {'ready':not errors and not pending,'errors':errors,'review_required':pending,'acknowledged':acknowledged,
            'target_sha256':th,'target_provenance':provenance,'phase_identity_assertions':assertions,'selection_sha256':digest(selection),
            'datasets_checked':len(inputs),'inputs':inputs,'weight_overrides':[{'path':r['path'],**r['weight_override']} for r in manifest.get('datasets',[]) if isinstance(r.get('weight_override'),dict)],
            'limitations':['Preflight checks input compatibility, not model accuracy or residual convergence.','ZPF uses driving-force tolerance; temperature/composition budgets are not automatically propagated.','Shared-source covariance is not estimated.'],'release_modified':False}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--selection',type=Path,required=True)
    a=ap.parse_args()
    try:r=check(a.bundle,a.selection)
    except Exception as exc:r={'ready':False,'errors':[{'code':'preflight_exception','message':type(exc).__name__+': '+str(exc)}]}
    print(json.dumps(r,indent=2));raise SystemExit(0 if r['ready'] else 1)
if __name__=='__main__':main()
