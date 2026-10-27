import copy
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import srsi_evidence_review as S
import dft_evidence_identity as I

@pytest.mark.parametrize('record', S.REVIEW['records'])
def test_classification_citation_and_protected_quantities(record):
    r=copy.deepcopy(record);before=copy.deepcopy(r);S.annotate([r])
    assert r['evidence_class']=='dft'
    assert r['cited_sources'][0]['doi']==S.NEW
    assert r['origin']['measured_in']=='doi:'+S.NEW
    allowed={'evidence_class','method','method_basis','origin','cited_sources','source_evidence_review','model_unlinked_reason'}
    assert {k:v for k,v in r.items() if k not in allowed}=={k:v for k,v in before.items() if k not in allowed}
    once=copy.deepcopy(r);S.annotate([r]);assert r==once

def test_changed_value_is_rejected():
    r=copy.deepcopy(S.REVIEW['records'][0]);r['value']=0
    with pytest.raises(ValueError,match='source evidence changed'):S.annotate([r])

def test_unreviewed_record_unchanged():
    r=copy.deepcopy(S.REVIEW['records'][0]);r['id']='UNREVIEWED';b=copy.deepcopy(r)
    S.annotate([r]);assert r==b

@pytest.mark.parametrize('pair', I.REVIEW['additional_reviews'][4]['pairs'])
def test_reclassified_repeated_results_do_not_gain_fitting_weight(pair):
    alias=copy.deepcopy(next(r for r in S.REVIEW['records'] if r['id']==pair['alias_id']))
    e=pair['record_evidence'][pair['canonical_id']]
    canonical=dict(id=pair['canonical_id'],paper=e['paper'],phases=[{'label':e['phase_label']}],temperature_K=e['temperature_K'],reference_state=e['reference_state'],unit='kJ/mol-atom',value=pair['value_kJ_mol_atom'],property='HM_FORM',evidence_class='dft')
    rs=[canonical,alias];S.annotate(rs);I.annotate(rs)
    assert [I.canonical(r) for r in rs]==[True,False]
    assert rs[0]['value']==rs[1]['value']==pair['value_kJ_mol_atom']
