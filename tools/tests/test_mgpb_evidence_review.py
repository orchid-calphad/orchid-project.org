import copy,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import mgpb_evidence_review as M
import uncertainty_budget as U
import build_espei_bundle as E

def records():
 p=M.REVIEW['pairs'][0]
 return [dict(id=rid,paper='doi:'+M.DOI,property='HM_FORM',unit='kJ/mol-atom',evidence_class='dft',temperature_K=0,temperature_printed='0 K',temperature_basis='printed',phases=[{'label':'Mg2Pb'}],flags={'extraction_accepted':True},uncertainty={'basis':'printed','value':.1},**v) for rid,v in p['records'].items()]

def test_preserves_numbers_and_separates_unit_identity():
 rs=records();old=copy.deepcopy(rs);M.annotate(rs)
 for a,b in zip(old,rs):
  for k in ['value','value_printed','unit_printed','unit','uncertainty','temperature_printed']:assert a[k]==b[k]
  assert b['evidence_class']=='experimental' and b['temperature_K'] is None
  assert b['temperature_K_stored']==0
 assert len({r['evidence_identity']['canonical_id'] for r in rs})==1
 snapshot=copy.deepcopy(rs);M.annotate(rs);assert snapshot==rs

def test_unresolved_conditions_excluded_from_calibration_and_fitting():
 rs=M.annotate(records())
 for r in rs:
  assert not U.usable(r)
  assert not E.lit_dft_candidate(r)
  assert E.gate_blocks(r)=='not_espei_ready:literature_gate'

def test_source_change_fails_closed():
 rs=records();rs[0]['value']=-99
 with pytest.raises(ValueError,match='source evidence changed'):M.annotate(rs)

def test_partial_pair_fails_closed():
 with pytest.raises(ValueError,match='Incomplete'):M.annotate(records()[:1])

def test_unrelated_rows_unchanged():
 rs=records()
 for i,r in enumerate(rs):r['id']='unrelated'+str(i)
 old=copy.deepcopy(rs);M.annotate(rs);assert rs==old
