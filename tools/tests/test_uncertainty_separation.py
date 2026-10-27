import copy
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import uncertainty_budget as U
import build_espei_bundle as E


def record(estimate, reported=None):
    return {'id':'ORC-TEST','uncertainty':estimate,'uncertainty_stored':reported,
            'unit':'kJ/mol-atom','uncertainty_budget':{'quantities':{'value':{'sigma_effective':5.63}}}}


def test_reported_and_scatter_are_separate_and_idempotent():
    r=record({'basis':'scatter_between_sources','value':5.63,'n_sources':2,'note':'ORC-OTHER'},
             {'basis':'printed','value':.24})
    expected_scatter=U.scatter_component(r,'value')
    U.separate_reported_uncertainty(r)
    assert r['uncertainty']=={'basis':'printed','value':.24}
    assert r['uncertainty_budget']['quantities']['value']['sigma_effective']==5.63
    assert U.scatter_component(r,'value')==expected_scatter
    assert E.printed_sigma(r,-12000)[0]==240
    old=copy.deepcopy(r);U.separate_reported_uncertainty(r);assert old==r


def test_composition_estimate_retained_and_missing_printed_not_invented():
    r=record({'basis':'printed_with_composition_offset','value':.788,'printed':.077,
              'sigma_composition_offset':.784,'note':'|dx| 0.0152 x slope 51.76'})
    expected=U.composition_component(r,'value');U.separate_reported_uncertainty(r)
    assert r['uncertainty']['value']==.077 and U.composition_component(r,'value')==expected
    r=record({'basis':'scatter_between_sources','value':8})
    U.separate_reported_uncertainty(r)
    assert r['uncertainty'] is None and E.printed_sigma(r,10)==(None,None)


def test_fitting_sigma_never_calls_legacy_estimate_printed():
    r=record({'basis':'scatter_between_sources','value':5.63},{'basis':'printed','value':.24})
    sigma,info=E.printed_sigma(r,-12000)
    assert sigma==240 and info['coverage']=='1sigma_assumed'
    r=record({'basis':'scatter_between_sources','value':5.63})
    assert E.printed_sigma(r,-12000)==(None,None)
