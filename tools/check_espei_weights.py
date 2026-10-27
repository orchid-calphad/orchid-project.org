#!/usr/bin/env python3
"""Check an ORCHID fitting bundle before applying any user weight overrides.
Usage: python check_weights.py --bundle /path/to/espei
Uses only the Python standard library; reads files without changing them.
"""
import argparse
import json
import math
from pathlib import Path


def check(bundle):
    errors = []
    index = json.loads((bundle / 'index.json').read_text())
    rows = index['datasets']; indexed = {r['path'] for r in rows}
    seen = set(); datasets = values = 0
    for row in rows:
        path = bundle / row['path']
        if not path.is_file():
            errors.append('missing indexed dataset: '+row['path']);continue
        d = json.loads(path.read_text());datasets += 1
        if d.get('output') == 'ZPF':continue
        u = d.get('orchid_uncertainty') or {}
        uv = u.get('values') or []
        if not uv:
            errors.append('missing uncertainty metadata: '+row['path']);continue
        default = (u.get('default_sigma') or {}).get('value')
        w = d.get('weight',1.0)
        weights = w if isinstance(w,list) else [w]*len(uv)
        if len(weights)!=len(uv):
            errors.append('weights/value count mismatch: '+row['path']);continue
        observations = (d.get('orchid') or {}).get('observations') or []
        for ob in observations:
            key=(row['folder'],ob['id'])
            if key in seen:errors.append('duplicate observation: '+str(key))
            seen.add(key)
        for weight,v in zip(weights,uv):
            sigma=v.get('sigma_effective_used')
            if sigma is None:sigma=(v.get('budget') or {}).get('sigma_effective')
            if sigma is None and v.get('sigma_basis')=='dft_calibrated':sigma=v.get('sigma')
            if not isinstance(sigma,(int,float)) or sigma<=0 or not math.isfinite(sigma):
                errors.append('missing/invalid effective sigma: '+row['path']);continue
            expected=round(default/sigma,6)
            if not math.isclose(weight,expected,rel_tol=0,abs_tol=0.00000051):
                errors.append('weight does not follow individual sigma: '+row['path'])
            if v.get('weight_used')!=weight:errors.append('weight metadata mismatch: '+row['path'])
            values+=1
    manifest=bundle/'UNCERTAINTY_SPLITS.json'
    if not manifest.exists():errors.append('missing split manifest')
    else:
        for split in json.loads(manifest.read_text())['splits']:
            if (bundle/split['parent_dataset']).exists():errors.append('stale parent: '+split['parent_dataset'])
            for child in split['children']:
                if child not in indexed:errors.append('unindexed child: '+child)
    for path in bundle.rglob('*.json'):
        rel=path.relative_to(bundle).as_posix()
        if rel in indexed:continue
        d=json.loads(path.read_text())
        if isinstance(d,dict) and 'output' in d:errors.append('unindexed fitting dataset: '+rel)
    return {'datasets_checked':datasets,'weighted_values_checked':values,'errors':errors}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parent)
    args=ap.parse_args();report=check(args.bundle)
    print(json.dumps(report,indent=2));raise SystemExit(bool(report['errors']))
