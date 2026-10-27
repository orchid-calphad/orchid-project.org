"""Regression protection for the Y5 known-issue patches (ORCHID-KI-20261027-v1), not a general rule:
(a) Ga-Sn Table 2 continuation rows extracted twice -> the page-3 copy gets duplicate_of the page-4 record;
(b) Ga-Sn metastable rows -> page-check ledger flag metastable_state (not_espei_ready), no value change;
(c) Fe-Y Table 3 Opt values at 973 K and Ag-Se Table 1 standard entropies at 298.15 K -> temperature from the page."""
import copy, json, sys, unittest
from pathlib import Path
TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
import human_review_corrections as h
import literature_gates as G
MAIN = json.loads((TOOLS / 'known_issues_corrections_20261027.json').read_text())
EXT_PATH = TOOLS / 'known_issues_corrections_20261027_extension.json'
EXT = json.loads(EXT_PATH.read_text()) if EXT_PATH.exists() else []
ALL = MAIN + EXT
KW = dict(evidence_id='ORCHID-KI-20261027-v1', gate='KI20261027', rule='source_adjudicated_known_issue',
          gap='known_issue_correction_requires_downstream_revalidation')
METASTABLE = ['ORC-NM3HGTM2DA', 'ORC-VLK3JFKMMT', 'ORC-QMU2IVHBUZ', 'ORC-ZA2QGJGERR', 'ORC-QDYBIJZXAK', 'ORC-436YUDCPRH',
              'ORC-MVGWPTTQLX', 'ORC-XKMT5GPSE7', 'ORC-6K3QBVA4HW', 'ORC-EBNBHQXQE3', 'ORC-IW5KWFEW7O', 'ORC-UJRYVPUOMY',
              'ORC-F6SDANDJJY', 'ORC-YZ73VJXLEF.3b206c', 'ORC-YZ73VJXLEF.c63806', 'ORC-YZ73VJXLEF.d4a478']
PRESERVED = ['value', 'unit', 'value_stored', 'uncertainty', 'cited_sources', 'pressure_Pa', 'overall_composition',
             'phases', 'phase_compositions', 'property', 'paper', 'pid', 'locator', 'evidence_class', 'merged_ids']


def fixtures(patches=ALL):
    rows = [dict(copy.deepcopy(p['before']), id=p['record_id'], flags={'espei_ready': False, 'extraction_accepted': True})
            for p in patches]
    for p in patches:   # the kept copy of each duplicate pair: a page-4 row that is itself not patched
        if 'duplicate_of' in p['after']:
            kept = dict(copy.deepcopy(p['before']), id=p['after']['duplicate_of'], flags={'extraction_accepted': True})
            kept['locator'] = dict(kept['locator'], page=4, row=p['evidence']['printed_row_page4'])
            rows.append(kept)
    return rows


def apply(rows, patches=None):
    if patches is None:
        return [h.apply(rows, MAIN, **KW), h.apply(rows, EXT, **KW)]
    return h.apply(rows, patches, **KW)


class KnownIssues(unittest.TestCase):
    def test_schema_and_scope(self):
        self.assertEqual(sum('duplicate_of' in p['after'] for p in MAIN), 31)
        self.assertEqual(sum(p['after'].get('temperature_K') == 973.0 for p in MAIN), 13)
        self.assertEqual(sum(p['after'].get('temperature_K') == 298.15 for p in MAIN), 16)
        for p in ALL:
            self.assertNotIn('flags', p['before'])          # gate-derived, rewritten later in the export
            self.assertEqual(p['evidence']['report_id'], 'ORCHID-KI-20261027-v1')
            self.assertEqual(p['visual_check']['verdict'], 'confirmed')
            self.assertTrue(set(p['after']) <= {'duplicate_of', 'temperature_K', 'temperature_basis'}, p['record_id'])
            for k in p['after']:
                self.assertIn(k, p['before'], 'every changed field is guarded: ' + k)
        self.assertEqual(len({p['record_id'] for p in ALL}), len(ALL))

    def test_duplicates_point_to_unpatched_page4_rows(self):
        dups = {p['record_id']: p for p in ALL if 'duplicate_of' in p['after']}
        for rid, p in dups.items():
            self.assertNotIn(p['after']['duplicate_of'], dups)   # no chains, no mutual pairs
            self.assertEqual(p['before']['locator']['page'], 3)   # the misplaced continuation copy
            self.assertIsNone(p['before']['duplicate_of'])
            self.assertTrue(p['evidence']['kept_pid'].startswith('doi:10.1007/bf02667485/p4/tab-2/'))
        self.assertEqual(len({p['after']['duplicate_of'] for p in dups.values()}), len(dups))

    def test_temperatures_only_where_printed(self):
        for p in ALL:
            t, b = p['after'].get('temperature_K'), p['before']
            if t is None:
                continue
            self.assertIsNone(b['temperature_K'])
            if t == 973.0:
                self.assertEqual(b['property_printed'], 'ΔHf,973'); self.assertEqual(p['after']['temperature_basis'], 'printed')
            elif t == 298.15:
                self.assertEqual(b['property'], 'SM'); self.assertEqual(b['paper'], 'doi:10.1007/s11669-018-0683-7')
            elif t == 0.0:
                self.assertEqual(b['property_printed'], 'ΔHf,0'); self.assertEqual(b['evidence_class'], 'dft')
            else:
                self.fail('unexpected temperature %r' % t)

    def test_apply_changes_only_patched_fields(self):
        rows = fixtures(); before = copy.deepcopy(rows); apply(rows)
        byid = {r['id']: r for r in rows}
        for p in ALL:
            r, old = byid[p['record_id']], next(x for x in before if x['id'] == p['record_id'])
            for k, v in p['after'].items():
                self.assertEqual(r[k], v)
            for k in PRESERVED:
                self.assertEqual(old.get(k), r.get(k), k)
            self.assertFalse(r['flags']['espei_ready'])
        kept = {p['after']['duplicate_of'] for p in ALL if 'duplicate_of' in p['after']}
        for old in before:
            if old['id'] in kept:
                self.assertEqual(byid[old['id']], old)   # the kept record is untouched
        self.assertEqual(len(rows), len(before))

    def test_replay_idempotent(self):
        rows = fixtures(); apply(rows); once = copy.deepcopy(rows); reports = apply(rows)
        self.assertEqual(rows, once)
        self.assertEqual([r['applied'] for r in reports], [[], []])

    def test_changed_guard_aborts_whole_batch(self):
        rows = fixtures(MAIN)
        fey = next(r for r in rows if r.get('property_printed') == 'ΔHf,973')
        fey['temperature_K'] = 298.15
        before = copy.deepcopy(rows)
        with self.assertRaises(ValueError):
            apply(rows, MAIN)
        self.assertEqual(rows, before)

    def test_missing_record_aborts_without_creating(self):
        rows = [r for r in fixtures(MAIN) if r['id'] != MAIN[0]['record_id']]
        before = copy.deepcopy(rows)
        with self.assertRaises(ValueError):
            apply(rows, MAIN)
        self.assertEqual(rows, before)

    def test_metastable_flags_in_ledger(self):
        visual = G.load_visual()
        for rid in METASTABLE:
            chk = visual.get('flag:' + rid)
            self.assertIsNotNone(chk, rid)
            self.assertIn({'code': 'metastable_state', 'effect': 'not_espei_ready'},
                          [{k: f[k] for k in ('code', 'effect')} for f in chk['flags']])
            self.assertTrue(chk.get('image_sha256')); self.assertEqual(chk.get('verdict'), 'confirmed')

    def test_metastable_flag_demotes_without_value_change(self):
        r = {'id': METASTABLE[0], 'value': None, 'temperature_K': 280.85, 'overall_composition': {'GA': 0.94, 'SN': 0.06},
             'flags': {'espei_ready': True, 'espei_certified': True, 'extraction_accepted': True}}
        old = copy.deepcopy(r)
        G.fail(r, 'page', 'metastable_state', 'not_espei_ready', 'metastable branch')
        self.assertFalse(r['flags']['espei_ready']); self.assertTrue(r['flags']['extraction_accepted'])
        for k in ('value', 'temperature_K', 'overall_composition'):
            self.assertEqual(r[k], old[k])

    def test_exporter_applies_the_batch(self):
        src = (TOOLS / 'export_literature_data.py').read_text()
        self.assertIn('known_issues_corrections_20261027.json', src)
        self.assertIn('ORCHID-KI-20261027-v1', src)


if __name__ == '__main__':
    unittest.main()
