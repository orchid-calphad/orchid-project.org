import copy,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import dft_evidence_identity as I

def records():
 p=I.REVIEW['pairs'][0]
 return [{'id':rid,'paper':'doi:'+I.REVIEW['source_doi'],'unit':'kJ/mol-atom','value':p['value_kJ_mol_atom'],'property':'HM_FORM','evidence_class':'dft','phases':[{'label':p['compound']}],'uncertainty':{'basis':'not printed'}} for rid in [p['canonical_id'],p['alias_id']]]
def test_values_preserved_and_one_identity_counted():
 r=records();old=copy.deepcopy(r);I.annotate(r)
 assert [I.canonical(x) for x in r]==[True,False]
 for a,b in zip(old,r):assert a=={k:v for k,v in b.items() if k!='evidence_identity'}
 snapshot=copy.deepcopy(r);I.annotate(r);assert snapshot==r

def test_changed_source_evidence_rejected():
 r=records();r[1]['value']=-1
 with pytest.raises(ValueError,match='evidence changed'):I.annotate(r)

def test_partial_pair_rejected():
 with pytest.raises(ValueError,match='Incomplete'):I.annotate(records()[:1])

def test_unreviewed_equal_values_not_merged():
 r=records()
 for i,x in enumerate(r):x['id']='UNREVIEWED-'+str(i)
 before=copy.deepcopy(r);I.annotate(r);assert r==before;assert all(I.canonical(x) for x in r)

def reporting_records(pair):
 return [dict(id=rid,paper=e['paper'],phases=[{'label':e['phase_label']}],temperature_K=e['temperature_K'],reference_state=copy.deepcopy(e['reference_state']),unit='kJ/mol-atom',value=pair['value_kJ_mol_atom'],property='HM_FORM',evidence_class='dft',uncertainty={'basis':'not printed'}) for rid,e in pair['record_evidence'].items()]

@pytest.mark.parametrize('pair', I.REVIEW['additional_reviews'][0]['pairs'])
def test_cross_paper_reporting_counts_once_without_changing_source_fields(pair):
 rs=reporting_records(pair);before=copy.deepcopy(rs);I.annotate(rs)
 assert sum(I.canonical(r) for r in rs)==1
 for a,b in zip(before,rs):assert a=={k:v for k,v in b.items() if k!='evidence_identity'}
 assert all(r['evidence_identity']['source_doi']=='10.1016/j.tca.2013.05.024' for r in rs)
 assert all('source_link_discrepancy' in r['evidence_identity'] for r in rs)
 snapshot=copy.deepcopy(rs);I.annotate(rs);assert snapshot==rs

def test_cross_paper_mismatched_state_rejected():
 rs=reporting_records(I.REVIEW['additional_reviews'][0]['pairs'][0]);rs[1]['temperature_K']=1000
 with pytest.raises(ValueError,match='evidence changed'):I.annotate(rs)

@pytest.mark.parametrize('pair', I.REVIEW['additional_reviews'][1]['pairs'])
def test_auce_original_and_later_report_preserve_distinct_conditions(pair):
 rs=reporting_records(pair);before=copy.deepcopy(rs);I.annotate(rs)
 assert [I.canonical(r) for r in rs]==[True,False]
 assert [r['temperature_K'] for r in rs]==[0.0,298.15]
 for a,b in zip(before,rs):assert a=={k:v for k,v in b.items() if k!='evidence_identity'}
 assert all('source_context_discrepancy' in r['evidence_identity'] for r in rs)
 assert all(r['evidence_identity']['source_doi']=='10.1016/j.calphad.2013.07.013' for r in rs)
 snap=copy.deepcopy(rs);I.annotate(rs);assert snap==rs

def test_mgsi_reports_share_primary_source_without_inventing_conditions():
 pair=I.REVIEW['additional_reviews'][2]['pairs'][0]
 rs=reporting_records(pair);rs[0].pop('temperature_K');before=copy.deepcopy(rs)
 I.annotate(rs)
 assert [I.canonical(r) for r in rs]==[True,False]
 assert 'temperature_K' not in rs[0] and rs[1]['temperature_K']==298.15
 for a,b in zip(before,rs):assert a=={k:v for k,v in b.items() if k!='evidence_identity'}
 assert rs[0]['evidence_identity']['source_doi']=='10.1016/j.intermet.2009.03.017'
 assert 'doi:'+rs[0]['evidence_identity']['canonical_reporting_source']['doi']==rs[0]['paper']
 snap=copy.deepcopy(rs);I.annotate(rs);assert snap==rs

def test_srsi_rounding_is_documented_without_replacing_reported_values():
 pair=I.REVIEW['additional_reviews'][3]['pairs'][0]
 rs=reporting_records(pair);before=copy.deepcopy(rs);I.annotate(rs)
 assert [I.canonical(r) for r in rs]==[True,False]
 for a,b in zip(before,rs):assert a=={k:v for k,v in b.items() if k!='evidence_identity'}
 for r in rs:
  assert r['value']==-46.9
  assert r['evidence_identity']['source_rounding']['original_source_value_kJ_mol_atom']==-46.89
  assert r['evidence_identity']['source_doi']=='10.1016/j.jallcom.2009.05.083'
  assert 'source_link_discrepancy' in r['evidence_identity']
 snap=copy.deepcopy(rs);I.annotate(rs);assert snap==rs
