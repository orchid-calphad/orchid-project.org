"""Regression protection for the source-verified Pd3Sn patch, not a general fill rule."""
import copy,json,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(TOOLS))
import human_review_corrections as h
PATCHES=json.loads((TOOLS/'table_context_corrections_20261027.json').read_text())
def fixtures():
 rows=[dict(copy.deepcopy(p['before']),id=p['record_id']) for p in PATCHES]
 # A following compound must not receive the preceding group's metadata.
 nextrow=copy.deepcopy(rows[-1]);nextrow['id']='ORC-NEXT-GROUP-TEST';nextrow['locator']['row']=7;nextrow['phases']=[{'label':'Pd2Sn'}]
 rows.append(nextrow);return rows
def apply(rows,patches=PATCHES):
 return h.apply(rows,patches,evidence_id='ORCHID-TC-20261027-v1',gate='TC20261027',rule='source_adjudicated_table_context')
class TableContext(unittest.TestCase):
 def test_first_and_continuation_rows_next_group_unchanged(self):
  rows=fixtures();before=copy.deepcopy(rows);apply(rows)
  self.assertEqual([r['locator']['row'] for r in rows[:-1]],[1,2,3,4,6])
  for old,r in zip(before[:-1],rows[:-1]):
   self.assertEqual(r['phases'][0]['label'],'Pd3Sn');self.assertEqual(r['structure']['space_group'],'Pm-3m')
   for field in ['value','unit','value_stored','uncertainty','cited_sources','temperature_K','pressure_Pa','overall_composition','merged_ids','merged_pids']:
    self.assertEqual(old.get(field),r.get(field),field)
  self.assertEqual(rows[-1],before[-1]);self.assertEqual(len(rows),len(before))
 def test_replay_idempotent(self):
  rows=fixtures();apply(rows);once=copy.deepcopy(rows);report=apply(rows)
  self.assertEqual(rows,once);self.assertEqual(report['applied'],[]);self.assertEqual(len(report['already_applied']),5)
 def test_conflicting_group_assignment_aborts_entire_batch(self):
  rows=fixtures();rows[2]['phases']=[{'label':'Pd2Sn'}];before=copy.deepcopy(rows)
  with self.assertRaises(ValueError):apply(rows)
  self.assertEqual(rows,before)
 def test_changed_row_reference_aborts(self):
  rows=fixtures();rows[1]['cited_sources']=[{'label':'different'}];before=copy.deepcopy(rows)
  with self.assertRaises(ValueError):apply(rows)
  self.assertEqual(rows,before)
 def test_missing_record_does_not_create_one(self):
  rows=fixtures()[1:];before=copy.deepcopy(rows)
  with self.assertRaises(ValueError):apply(rows)
  self.assertEqual(rows,before)
if __name__=='__main__':unittest.main()
