"""Source-checked Sr-Si classification and citation repair, preserving reported quantities."""
import copy
import json
from pathlib import Path

REVIEW = json.loads(Path(__file__).with_suffix('.json').read_text())
OLD = REVIEW['existing_resolved_doi']
NEW = REVIEW['verified_original_doi']


def annotate(records):
    expected = {r['id']: r for r in REVIEW['records']}
    for r in records:
        if r['id'] not in expected:
            continue
        e = expected[r['id']]
        for key in ('paper', 'property', 'value', 'unit', 'phases', 'temperature_K', 'reference_state', 'uncertainty'):
            if r.get(key) != e.get(key):
                raise ValueError('Sr-Si source evidence changed: '+r['id']+' '+key)
        done = r.get('source_evidence_review', {}).get('review_id') == REVIEW['review_id']
        method = 'first-principles '+REVIEW['methods'][r['id']]
        if done:
            if r.get('evidence_class') != 'dft' or r.get('method') != method:
                raise ValueError('Sr-Si reviewed classification changed: '+r['id'])
        else:
            if r.get('evidence_class') != e['evidence_class']:
                raise ValueError('Unexpected Sr-Si classification: '+r['id'])
            before = {k:copy.deepcopy(r.get(k)) for k in ('evidence_class','origin','method','method_basis','model_unlinked_reason')}
            r['evidence_class'] = 'dft'
            r['method'] = method
            r['method_basis'] = 'GGA/LDA printed in Table 2, PDF page 2'
            r.pop('model_unlinked_reason', None)
            r['source_evidence_review'] = {'review_id':REVIEW['review_id'], 'source_doi':REVIEW['reporting_doi'],
                'pdf_page':2, 'table':2, 'before':before,
                'fitting_status':'repeated_report_use_canonical_evidence_identity',
                'basis':REVIEW['basis']}
        cites = [c for c in r.get('cited_sources', []) if str(c.get('label')) == '17']
        if len(cites) != 1 or cites[0].get('reference') != e['cited_sources'][0]['reference']:
            raise ValueError('Sr-Si printed reference changed: '+r['id'])
        c = cites[0]
        if c.get('doi') not in (OLD, NEW) or c.get('paper') != 'doi:'+c['doi']:
            raise ValueError('Unexpected Sr-Si source identity: '+r['id'])
        c.update(doi=NEW, paper='doi:'+NEW, kind='calculation', kind_original=e['cited_sources'][0]['kind'],
                 doi_original=OLD, doi_status='corrected_from_printed_reference', doi_correction=REVIEW['review_id'])
        origin = r.setdefault('origin', {})
        if origin.get('measured_in') not in ('doi:'+OLD, 'doi:'+NEW):
            raise ValueError('Unexpected Sr-Si original source: '+r['id'])
        origin.update(measured_in='doi:'+NEW, measured_in_original='doi:'+OLD,
                      source_reference_correction=REVIEW['review_id'],
                      **{'class':'computed','basis':'printed_GGA_LDA_cited_calculation','own_work_of_paper':False})
    return records


def reconcile_papers(records, papers):
    pk = 'doi:'+NEW
    papers.setdefault(pk, {'id':pk, 'doi':NEW, 'doi_url':'https://doi.org/'+NEW,
        'title':'Thermodynamics, lattice stability and defect structure of strontium silicides via first-principles calculations',
        'authors':['A. Garay','M.E. Williams','G. Trápaga','R. Arróyave'], 'year':2009,
        'container':'Journal of Alloys and Compounds', 'kind':'cited_work',
        'citation':'A. Garay, M.E. Williams, G. Trápaga, R. Arróyave, Thermodynamics, lattice stability and defect structure of strontium silicides via first-principles calculations, Journal of Alloys and Compounds 484, 822–831 (2009).',
        'metadata_source':'verified_source_pdf_and_printed_reference',
        'roles':[], 'source_kinds':[]})
    for doi in (OLD, NEW):
        pk = 'doi:'+doi
        p = papers[pk]
        reporting = [r for r in records if r.get('paper') == pk]
        cited = [r for r in records if any(c.get('paper') == pk for c in r.get('cited_sources', []))]
        linked = [r for r in records if r.get('origin', {}).get('measured_in') == pk]
        p['systems'] = sorted({r['system'] for r in reporting+cited+linked})
        roles = set(p.get('roles', []))-{'cited_source'}
        kinds = set(p.get('source_kinds', []))-{'cited_work'}
        if cited:
            roles.add('cited_source');kinds.add('cited_work')
        p['roles'],p['source_kinds'] = sorted(roles),sorted(kinds)
    return papers
