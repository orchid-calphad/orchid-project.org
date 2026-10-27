import copy
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import boride_evidence_review as B
import build_espei_bundle as E
import uncertainty_budget as U


def test_preservation_idempotence_and_fitting_exclusion():
    records=copy.deepcopy(B.REVIEW['records'])
    old=copy.deepcopy(records)
    B.annotate(records)
    once=copy.deepcopy(records)
    B.annotate(records)
    assert records==once
    for a,b in zip(old,records):
        allowed={'evidence_class','method','method_basis','origin','source_evidence_review','gate_failures'}
        assert {k:v for k,v in a.items() if k not in allowed}=={k:v for k,v in b.items() if k not in allowed}
        assert not E.lit_dft_candidate(b)
        assert not U.usable(b)
    assert [r['evidence_class'] for r in records]==['experimental','dft','dft']


@pytest.mark.parametrize('field,value',[('value',0),('uncertainty',{}),('temperature_K',1473),('reference_state',{}),('cited_sources',[])])
def test_guard_rejects_source_changes(field,value):
    records=copy.deepcopy(B.REVIEW['records']);records[0][field]=value
    with pytest.raises(ValueError):B.annotate(records)


def test_identity_exclusion_preserves_both_records_and_is_idempotent():
    records=copy.deepcopy(B.REVIEW['identity_review']['records']);old=copy.deepcopy(records)
    B.annotate(records);once=copy.deepcopy(records);B.annotate(records)
    assert records==once and len(records)==2
    for a,b in zip(old,records):
        assert {k:v for k,v in b.items() if k not in {'source_evidence_review','gate_failures'}}==a
        assert not E.lit_dft_candidate(b) and not U.usable(b)
        assert 'evidence_identity' not in b
        assert b['source_evidence_review']['identity_status']=='unresolved_not_asserted_duplicate'


@pytest.mark.parametrize('field,value',[('value',0),('cited_sources',[]),('uncertainty',{}),('temperature_K',0)])
def test_identity_guard_rejects_source_changes(field,value):
    records=copy.deepcopy(B.REVIEW['identity_review']['records']);records[0][field]=value
    with pytest.raises(ValueError):B.annotate(records)
