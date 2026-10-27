"""Guarded source-adjudicated corrections. No network access or publication.

Apply to fully exported records, then regenerate derived views and uncertainty
budgets. A conflict aborts the whole batch before any record is modified.
"""
import copy

EVIDENCE_ID = 'ORCHID-HR-20261026-v1'


def apply(records, patches, evidence_id=EVIDENCE_ID, gate='HR20261026', rule='source_adjudicated_human_review',
          gap='human_review_correction_requires_downstream_revalidation'):
    index = {r['id']: r for r in records}
    pending, already, errors = [], [], []
    for patch in patches:
        rec = index.get(patch['record_id'])
        if rec is None:
            errors.append(patch['record_id'] + ': missing record')
            continue
        if all(rec.get(k) == v for k, v in patch['after'].items()):
            already.append(patch['record_id'])
            continue
        differences = [k for k, v in patch['before'].items() if rec.get(k) != v]
        if differences:
            errors.append(patch['record_id'] + ': changed guarded fields ' + ', '.join(differences))
            continue
        pending.append((rec, patch))
    if errors:
        raise ValueError('Correction batch aborted without changes:\n' + '\n'.join(errors))
    for rec, patch in pending:
        for field, after in patch['after'].items():
            before = copy.deepcopy(rec.get(field))
            if before == after:
                continue
            rec[field] = copy.deepcopy(after)
            rec.setdefault('corrections', []).append({
                'gate': gate, 'field': field, 'before': before,
                'after': copy.deepcopy(after), 'rule': rule,
                'evidence': copy.deepcopy(patch['evidence']),
                'visual_check': copy.deepcopy(patch['visual_check'])})
        # Old model agreement cannot certify a corrected state or provenance.
        flags = rec.setdefault('flags', {})
        flags['espei_ready'] = False
        flags['espei_certified'] = False
        flags['state_fully_specified'] = False
        gaps = flags.setdefault('specification_gaps', [])
        if gap not in gaps:
            gaps.append(gap)
        if flags.get('plausibility') is not None:
            rec.setdefault('plausibility_before_correction', copy.deepcopy(flags['plausibility']))
            flags.pop('plausibility')
        flags['calphad_check'] = 'not_available'
        if patch.get('residual_concern'):
            rec.setdefault('gate_failures', []).append({
                'gate': gate, 'code': 'human_review_residual_concern',
                'effect': 'not_espei_ready', 'detail': patch['residual_concern']})
    return {'evidence_id': evidence_id, 'applied': [p['record_id'] for _, p in pending],
            'already_applied': already, 'fields_changed': sum(
                sum(p['before'].get(k) != v for k, v in p['after'].items()) for _, p in pending)}
