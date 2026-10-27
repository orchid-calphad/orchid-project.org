#!/usr/bin/env python3
"""Read-only audit of calculated diagrams. Reports never alter release status.

Boundary probes preserve composition-set multiplicity. Grid probes compare phase
sets because the stored grid labels do not distinguish repeated composition sets.
Passing sampled probes does not establish complete or experimentally correct maps.
"""
import argparse, collections, hashlib, json, math, multiprocessing as mp, os
from pathlib import Path
import signal, time, re
from compute_binary_diagrams import prepare_source, source_phase_selection
from diagram_temperature_extent import temperature_extent

def phase_name(label):return re.sub(r'#\d+$','',label)

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
def structural(p,tdbs):
    d=json.loads(p.read_text());issues=[];warnings=[];phases=set(d.get('phases',[]))
    r={'file':p.name,'diagram_sha256':digest(p),'source_sha256':d.get('sha256'),'issues':issues,'warnings':warnings,'probes':[]}
    source=tdbs/d.get('tdb','')
    if not source.is_file():issues.append('source_missing')
    elif digest(source)!=d.get('sha256'):issues.append('source_hash_mismatch')
    lo,hi=d.get('T_range',[float('nan')]*2)
    if not (math.isfinite(lo) and math.isfinite(hi) and 0<lo<hi):issues.append('invalid_temperature_range')
    def coordinate(x,t):return math.isfinite(x) and -1e-5<=x<=1+1e-5 and math.isfinite(t) and lo-.02<=t<=hi+.02
    documented_extent=False
    if 'temperature_extent' in d:
        try:
            documented_extent=d['temperature_extent']==temperature_extent(d)
        except (ValueError,TypeError,KeyError):
            pass
        if not documented_extent:issues.append('temperature_extent_metadata_mismatch')
        elif d['temperature_extent']['outside_requested_range']:
            warnings.append('documented_boundary_extent_outside_requested_range')
    def boundary_coordinate(x,t):
        return (math.isfinite(x) and -1e-5<=x<=1+1e-5 and math.isfinite(t) and t>0
                and (documented_extent or lo-.02<=t<=hi+.02))
    choices=[]
    invariant_temperatures=[inv['T'] for inv in d.get('invariants',[]) if math.isfinite(inv['T'])]
    for ri,reg in enumerate(d.get('regions',[])):
        for branch in reg:
            xs,ts=branch.get('x',[]),branch.get('T',[])
            if phase_name(branch.get('phase','')) not in phases:issues.append('unknown_boundary_phase')
            if len(xs)!=len(ts):issues.append('boundary_array_length_mismatch')
            if not all(boundary_coordinate(x,t) for x,t in zip(xs,ts)):issues.append('invalid_boundary_coordinate')
        if len(reg)!=2 or len(reg[0]['x'])!=len(reg[1]['x']):continue
        candidates=[]
        for i,(a,b) in enumerate(zip(reg[0]['x'],reg[1]['x'])):
            t=reg[0]['T'][i];x=(a+b)/2
            if (abs(t-reg[1]['T'][i])<.03 and abs(a-b)>.002 and .002<x<.998 and coordinate(x,t)
                and all(abs(t-it)>.1 for it in invariant_temperatures)):
                candidates.append({'kind':'boundary','region':ri,'T_K':t,'x':x,'expected':sorted(phase_name(z['phase']) for z in reg),'width':abs(a-b)})
        if candidates:choices.append(max(candidates,key=lambda c:c['width']))
    for inv in d.get('invariants',[]):
        for point in inv.get('points',[]):
            if phase_name(point.get('phase','')) not in phases:issues.append('unknown_invariant_phase')
            if not coordinate(point['x'],inv['T']):issues.append('invalid_invariant_coordinate')
    if not d.get('regions') or not any(b.get('x') for reg in d.get('regions',[]) for b in reg):warnings.append('empty_boundaries_requires_review')
    # Widest region first, then distinct phase assemblages, up to two probes.
    choices.sort(key=lambda c:c['width'],reverse=True);chosen=[]
    for c in choices:
        if not chosen or c['expected']!=chosen[0]['expected']:
            chosen.append(c)
        if len(chosen)==2:break
    if len(chosen)<2 and len(choices)>1:
        other=max(choices,key=lambda c:abs(c['T_K']-chosen[0]['T_K']))
        if other not in chosen:chosen.append(other)
    grid=d.get('grid')
    if grid:
        rows=grid.get('rows',[]);labels=grid.get('labels',[]);nx=grid.get('x_points',0)
        for label in labels:
            if label and any(ph not in phases for ph in label.split('+')):issues.append('unknown_grid_phase')
        for row in rows:
            if sum(run[1] for run in row)!=nx:issues.append('grid_row_length_mismatch')
            if any(not isinstance(run[1],int) or run[1]<=0 or not isinstance(run[0],int) or not 0<=run[0]<len(labels) for run in row):issues.append('invalid_grid_run')
        if rows and nx>2:
            for tf,xf in [(.25,.33),(.75,.67)]:
                ti=round((len(rows)-1)*tf);xi=round((nx-1)*xf);cursor=0;label=None
                for labelidx,n in rows[ti]:
                    if cursor<=xi<cursor+n and 0<=labelidx<len(labels):label=labels[labelidx];break
                    cursor+=n
                t=grid['T0']+ti*grid['dT'];x=xi*grid['x_step']
                if not coordinate(x,t):issues.append('invalid_grid_coordinate')
                elif label:chosen.append({'kind':'grid','T_K':t,'x':x,'expected':sorted(set(label.split('+')))})
                else:warnings.append('sampled_grid_cell_has_no_phase')
    else:
        warnings.append('no_equilibrium_grid')
        if len(choices)>2:
            mid=choices[len(choices)//2]
            if mid not in chosen:chosen.append(mid)
    if not chosen:warnings.append('no_checkable_probe')
    r['probes']=chosen;r['issues']=sorted(set(issues));r['warnings']=sorted(set(warnings))
    return r

class AuditTimeout(Exception):pass
def alarm(*a):raise AuditTimeout('diagram audit time limit')
def evaluate(task):
    path,tdbs,out,seconds,pdens=task;p=Path(path);dest=Path(out)/'equilibrium'/p.name
    t0=time.monotonic();r={'file':p.name,'checks':[],'issues':[],'warnings':[]}
    signal.signal(signal.SIGALRM,alarm);signal.alarm(seconds)
    try:
        import numpy as np
        from pycalphad import Database,equilibrium,variables as v
        r=structural(p,Path(tdbs));r['checks']=[]
        d=json.loads(p.read_text());text,_,adjustments=prepare_source((Path(tdbs)/d['tdb']).read_bytes())
        if d.get('parsed_text_sha256') and hashlib.sha256(text.encode()).hexdigest()!=d['parsed_text_sha256']:r['issues'].append('parsed_source_hash_mismatch')
        db=Database(text);active=sorted(db.phases)
        try:
            defaults,selection=source_phase_selection(text,active)
        except ValueError as exc:
            if d.get('phase_selection'):raise
            # Unused source defaults do not define the stored calculation domain.
            defaults,selection=active,{'excluded':[]}
            r['warnings'].append('unresolved_unused_source_defaults')
            r['source_defaults_error']=str(exc)
        if d.get('phase_selection'):
            active=defaults
            if selection!=d['phase_selection']:r['issues'].append('phase_selection_provenance_mismatch')
        elif selection['excluded']:r['warnings'].append('source_default_exclusions_not_applied')
        if set(active)!=set(d['phases']):r['issues'].append('active_phase_list_mismatch')
        if 'source_hash_mismatch' in r['issues'] or 'active_phase_list_mismatch' in r['issues']:raise ValueError('Cannot validate against mismatched source or phase domain')
        comps=d['components']+(['VA'] if 'VA' in db.elements else [])
        for probe in r['probes']:
            c=dict(probe)
            try:
                eq=equilibrium(db,comps,active,{v.X(d['x_component']):probe['x'],v.T:probe['T_K'],v.P:101325,v.N:1},calc_opts={'pdens':pdens})
                labels=eq.Phase.values.flatten();amounts=eq.NP.values.flatten()
                got=sorted(str(ph) for ph,amt in zip(labels,amounts) if str(ph) and np.isfinite(amt) and amt>1e-8)
                if probe['kind']=='grid':got=sorted(set(got))
                c.update(actual=got,finite_GM=bool(np.all(np.isfinite(eq.GM.values))))
                c['agrees']=c['finite_GM'] and got==probe['expected']
            except AuditTimeout:raise
            except Exception as exc:c.update(agrees=False,error=type(exc).__name__+': '+str(exc)[:350])
            r['checks'].append(c)
        if not r['checks']:r['result']='not_checkable'
        elif any(not c['agrees'] for c in r['checks']):r['result']='requires_review'
        elif r['issues']:r['result']='structural_issue'
        else:r['result']='sampled_checks_pass'
    except AuditTimeout:r['result']='timeout'
    except Exception as exc:r.update(result='engine_or_source_error',error=type(exc).__name__+': '+str(exc)[:500])
    finally:signal.alarm(0)
    r['seconds']=round(time.monotonic()-t0,2);r['protocol']={'pdens':pdens,'phase_amount_threshold':1e-8,'invariant_avoidance_K':.1,'composition_set_suffix_normalized':True,'timeout_seconds':seconds,'source_of_pressure':'mapping protocol, 101325 Pa','scope':'Sampled consistency only; no experimental or completeness certification.'}
    write(dest,r);return p.name,r['result'],r['seconds']

def main():
    ap=argparse.ArgumentParser();h=Path.home()
    ap.add_argument('--data',type=Path,default=h/'projects/orchid-data');ap.add_argument('--tdbs',type=Path,default=h/'projects/open-tdb/databases/sgte-binary');ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=6);ap.add_argument('--timeout',type=int,default=60);ap.add_argument('--pdens',type=int,default=1000);ap.add_argument('--structural-only',action='store_true');ap.add_argument('--resume',action='store_true');ap.add_argument('--limit',type=int)
    a=ap.parse_args();paths=[p for p in sorted((a.data/'phase-diagrams/sgte').glob('*.json')) if json.loads(p.read_text()).get('status')=='ok']
    if a.limit:paths=paths[:a.limit]
    structures=[]
    for p in paths:
        try:r=structural(p,a.tdbs)
        except Exception as exc:r={'file':p.name,'diagram_sha256':digest(p),'issues':['structural_audit_exception'],'error':repr(exc)}
        structures.append(r)
    write(a.out/'structural.json',{'diagrams':len(paths),'issue_counts':dict(collections.Counter(i for r in structures for i in r['issues'])),'warning_counts':dict(collections.Counter(i for r in structures for i in r.get('warnings',[]))),'records':structures})
    print('Structural audit:',len(paths),'diagrams,',sum(bool(r['issues']) for r in structures),'flagged',flush=True)
    if a.structural_only:return
    tasks=[]
    for p in paths:
        prev=a.out/'equilibrium'/p.name
        if a.resume and prev.exists():
            r=json.loads(prev.read_text())
            if r.get('diagram_sha256')==digest(p) and r.get('source_sha256')==digest(a.tdbs/json.loads(p.read_text())['tdb']) and r['protocol']['pdens']==a.pdens:continue
        tasks.append((str(p),str(a.tdbs),str(a.out),a.timeout,a.pdens))
    with mp.Pool(a.workers,maxtasksperchild=25) as pool:
        for i,res in enumerate(pool.imap_unordered(evaluate,tasks),1):
            if i%25==0 or i==len(tasks):print(i,'/',len(tasks),res,flush=True)
    rows=[json.loads((a.out/'equilibrium'/p.name).read_text()) for p in paths]
    write(a.out/'SUMMARY.json',{'diagrams':len(paths),'results':dict(collections.Counter(r['result'] for r in rows)),'checks':sum(len(r['checks']) for r in rows),'release_modified':False})
if __name__=='__main__':main()
