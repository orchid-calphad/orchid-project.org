"""Source-adjudicated placeholder recovery: guards, provenance, and gate safety."""
import copy
import json
from pathlib import Path
import sys
import unittest
TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
from human_review_corrections import apply
PATCHES = json.loads((TOOLS / 'placeholder_recovery_corrections_20261027.json').read_text())
KW = dict(evidence_id='ORCHID-PR-20261027-v1', gate='PR20261027',
          rule='source_adjudicated_placeholder_recovery',
          gap='placeholder_recovery_requires_downstream_revalidation')
def fixtures(patches):
    return [dict(copy.deepcopy(p['before']), id=p['record_id'],
                 flags={'extraction_accepted': False, 'espei_ready': True,
                        'espei_certified': True, 'plausibility': {'result': 'pass'}})
            for p in patches]
class PlaceholderRecovery(unittest.TestCase):
    def test_phase_only_and_source_evidence(self):
        self.assertEqual(len(PATCHES), 1267)
        self.assertEqual(len({p['record_id'] for p in PATCHES}), len(PATCHES))
        for p in PATCHES:
            self.assertEqual(set(p['after']), {'phases'})
            self.assertNotIn('flags', p['before'])
            self.assertTrue(p['after']['phases'])
            self.assertTrue(all(x['label_source']=='printed_source_page' for x in p['after']['phases']))
            self.assertLessEqual(len(p['evidence']['printed_phase_text']),500)
            self.assertEqual(len(p['evidence']['rendered_page_sha256']),64)
    def test_apply_replay_and_no_promotion(self):
        rows=fixtures(PATCHES);apply(rows,PATCHES,**KW)
        for row in rows:
            self.assertFalse(row['flags']['extraction_accepted'])
            self.assertFalse(row['flags']['espei_ready'])
            self.assertFalse(row['flags']['espei_certified'])
            self.assertNotIn('plausibility',row['flags'])
        saved=copy.deepcopy(rows);result=apply(rows,PATCHES,**KW)
        self.assertEqual(rows,saved)
        self.assertEqual(len(result['already_applied']),1267)
    def test_numeric_guard_aborts_all(self):
        rows=fixtures(PATCHES[:2]);rows[-1]['value']=987654321
        saved=copy.deepcopy(rows)
        with self.assertRaises(ValueError):apply(rows,PATCHES[:2],**KW)
        self.assertEqual(rows,saved)
    def test_value_guard_aborts_all(self):
        # cited_sources left the before-guard at the refresh (enriched later in the export); the printed value still guards
        rows=fixtures(PATCHES[:2]);rows[-1]['value']='altered'
        saved=copy.deepcopy(rows)
        with self.assertRaises(ValueError):apply(rows,PATCHES[:2],**KW)
        self.assertEqual(rows,saved)
    def test_missing_record_aborts_without_creation(self):
        rows=fixtures(PATCHES[:1]);saved=copy.deepcopy(rows)
        with self.assertRaises(ValueError):apply(rows,PATCHES[:2],**KW)
        self.assertEqual(rows,saved)
    def test_unlisted_neighbor_untouched(self):
        rows=fixtures(PATCHES[:1]);neighbor=copy.deepcopy(rows[0]);neighbor['id']='UNLISTED-NEIGHBOR'
        rows.append(neighbor);saved=copy.deepcopy(neighbor);apply(rows,PATCHES[:1],**KW)
        self.assertEqual(rows[1],saved)
    def test_continuation_and_withheld_rows(self):
        excluded={'ORC-SWBBBB3ADA','ORC-D5SSCLQRHX','ORC-MM4BTTUEOP','ORC-S4ZIORSOPU'}
        self.assertFalse(excluded.intersection(p['record_id'] for p in PATCHES))
        continuation=[p for p in PATCHES if p['evidence']['review_group']==1 and p['before']['locator']['row']>185]
        self.assertTrue(continuation)
        self.assertTrue(all(p['evidence']['pdf_page']==9 for p in continuation))
if __name__=='__main__':unittest.main()
