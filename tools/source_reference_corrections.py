"""Guarded, source-page-verified citation corrections; numerical fields are immutable."""
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
        if r.get('paper') != 'doi:' + REVIEW['reporting_doi'] or any(
            r.get(k) != e[k] for k in ('property', 'value', 'unit', 'evidence_class')
        ):
            raise ValueError('Source-reference evidence changed: ' + r['id'])
        citations = [c for c in r.get('cited_sources', []) if str(c.get('label')) == REVIEW['reference_label']]
        if len(citations) != 1:
            raise ValueError('Expected one reference 16: ' + r['id'])
        c = citations[0]
        if c.get('reference') != e['cited_sources'][0]['reference'] or c.get('doi') not in (OLD, NEW):
            raise ValueError('Printed source-reference changed: ' + r['id'])
        if c.get('paper') != 'doi:' + c['doi']:
            raise ValueError('Inconsistent citation identity: ' + r['id'])
        c.update(doi=NEW, paper='doi:' + NEW,
                 doi_status='corrected_from_printed_reference', doi_correction=REVIEW['review_id'],
                 doi_original=OLD)
        origin = r.get('origin', {})
        if origin.get('measured_in') not in ('doi:' + OLD, 'doi:' + NEW):
            raise ValueError('Unexpected original source: ' + r['id'])
        origin['measured_in'] = 'doi:' + NEW
        origin['measured_in_original'] = 'doi:' + OLD
        origin['source_reference_correction'] = REVIEW['review_id']
    return records


def reconcile_papers(records, papers):
    """Recalculate affected paper memberships without removing the valid La-B/La-Bi paper."""
    for doi in (OLD, NEW):
        pk = 'doi:' + doi
        p = papers[pk]
        reporting = [r for r in records if r.get('paper') == pk]
        cited = [r for r in records if any(c.get('paper') == pk for c in r.get('cited_sources', []))]
        linked = [r for r in records if r.get('origin', {}).get('measured_in') == pk]
        p['systems'] = sorted({r['system'] for r in reporting + cited + linked})
        roles = set(p.get('roles', [])) - {'cited_source'}
        kinds = set(p.get('source_kinds', [])) - {'cited_work'}
        if cited:
            roles.add('cited_source')
            kinds.add('cited_work')
        p['roles'], p['source_kinds'] = sorted(roles), sorted(kinds)
    return papers
