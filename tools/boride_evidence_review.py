"""Source-checked boride classifications; preserve reported quantities and conditions."""
import copy
import json
from pathlib import Path

REVIEW = json.loads(Path(__file__).with_suffix('.json').read_text())
GATE = 'boride_primary_source_state_review_pending'
IDENTITY_GATE = 'boride_source_identity_unresolved'


def fitting_blocked(record):
    return any(g.get('code') in {GATE, IDENTITY_GATE} for g in record.get('gate_failures', []))


def annotate(records):
    expected = {r['id']: r for r in REVIEW['records']}
    for r in records:
        if r['id'] not in expected:
            annotate_identity(r)
            continue
        old = expected[r['id']]
        decision = REVIEW['decisions'][r['id']]
        for key in ('paper', 'property', 'value', 'value_printed', 'unit', 'phases',
                    'temperature_K', 'reference_state', 'uncertainty', 'cited_sources'):
            if r.get(key) != old.get(key):
                raise ValueError('Boride source evidence changed: '+r['id']+' '+key)
        if r.get('source_evidence_review', {}).get('review_id') == REVIEW['review_id']:
            if (r.get('evidence_class') != decision['evidence_class'] or
                r.get('method') != decision['method'] or not fitting_blocked(r)):
                raise ValueError('Boride reviewed state changed: '+r['id'])
            continue
        if r.get('evidence_class') != old['evidence_class'] or r.get('method') != old['method']:
            raise ValueError('Unexpected boride classification: '+r['id'])
        before = {k:copy.deepcopy(r.get(k)) for k in ('evidence_class','method','method_basis','origin')}
        r.update(evidence_class=decision['evidence_class'],method=decision['method'],method_basis=decision['basis'])
        r['origin'] = dict(r.get('origin', {}), **{'class':'computed' if decision['evidence_class']=='dft' else 'experimental', 'basis':'source_checked_boride_classification'})
        r['source_evidence_review'] = {'review_id':REVIEW['review_id'],'before':before,
            'source_doi':REVIEW['reporting_doi'],'pdf_page':25,'table':2,'row':r['locator']['row'],
            'decision':copy.deepcopy(decision),'fitting_status':'excluded_pending_primary_source_state_verification'}
        r.setdefault('gate_failures', []).append({'code':GATE,'effect':'not_espei_ready',
            'reason':'Evidence class corrected; preserve printed conditions and uncertainties. Complete original-source numerical/reference-state verification before fitting or calibration.'})
    return records


def reconcile_papers(records, papers):
    # Printed citation labels remain intact; verified source links live in the review.
    return papers


def annotate_identity(r):
    review = REVIEW['identity_review']
    expected = {x['id']: x for x in review['records']}
    if r['id'] not in expected:
        return
    old = expected[r['id']]
    # Only the review and gate may change; all source fields remain immutable.
    allowed = {'source_evidence_review', 'gate_failures', 'uncertainty_budget', 'fold_lso', 'fold_lpo', 'chemical_system', 'paper_id'}   # the folds are added after this annotation
    if {k:v for k,v in r.items() if k not in allowed} != {k:v for k,v in old.items() if k not in allowed}:
        raise ValueError('LuB12 source evidence changed: '+r['id'])
    if r.get('source_evidence_review', {}).get('review_id') == review['review_id']:
        if not any(g.get('code') == IDENTITY_GATE for g in r.get('gate_failures', [])):
            raise ValueError('LuB12 exclusion gate missing: '+r['id'])
        return
    if r.get('source_evidence_review'):
        raise ValueError('Unexpected prior LuB12 review: '+r['id'])
    r['source_evidence_review'] = {'review_id':review['review_id'],
        'source_doi':review['reporting_doi'], 'pdf_pages':review['source_pages'],
        'table':2, 'row':r['locator']['row'], 'basis':review['basis'],
        'fitting_status':'excluded_pending_independent_source_identity',
        'identity_status':'unresolved_not_asserted_duplicate',
        'before':{'gate_failures':copy.deepcopy(r.get('gate_failures'))},
        'readmission':review['readmission']}
    r.setdefault('gate_failures', []).append({'code':IDENTITY_GATE, 'effect':'not_espei_ready',
        'reason':review['decision']})
