"""Regression coverage for unequal uncertainty groups, point alignment and source preservation."""
import copy
from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_espei_bundle as E


def fixture(sigmas=(0.5, 1.5), solver=False):
    pub = {f'ORC-{i}': {'id': f'ORC-{i}', 'unit': 'kJ/mol-atom', 'paper': 'doi:shared',
                       'uncertainty': {'basis': 'printed', 'value': s},
                       'uncertainty_budget': {'quantities': {'value': {'sigma_effective': s}}}}
           for i, s in enumerate(sigmas)}
    recs = [(None, r) for r in pub.values()]
    d = {'output': 'HM_FORM' if solver else 'HMR', 'values': [[[-1000., -2000.]]],
         'conditions': {'P': 101325, 'T': 500., 'X_B': [.1, .9]}, 'tags': [],
         'comment': 'ORCHID observation ids (in value order): ORC-0 ORC-1',
         'orchid': {'observations': [{'id': r['id']} for _, r in recs], 'source_dataset': 'source.json'}}
    if solver:
        d['solver'] = {'sublattice_configurations': [[['A']], [['B']]],
                       'sublattice_occupancies': [[[1.]], [[1.]]]}
        d['conditions'].pop('X_B')
    weight, block, tags, stat = E.uncertainty_block(d, recs, 'experimental', False, 500.)
    d.update(weight=weight, orchid_uncertainty=block)
    rel = 'experimental/A-B/Test.hmr.json'
    c = {'datasets': {rel: d}, 'rows': {rel: {'evidence_class': 'experimental'}},
         '_recs': {rel: list(pub.values())}, '_ustat': {rel: stat}}
    return c, pub, rel


def test_split_preserves_values_conditions_budgets_and_source():
    c, pub, rel = fixture(); original = copy.deepcopy(c['datasets'][rel])
    splits = E.split_unequal_uncertainties(c, pub)
    assert len(splits) == 1 and len(c['datasets']) == 2 and rel not in c['datasets']
    for i, child in enumerate(c['datasets'].values()):
        assert child['values'] == [[[original['values'][0][0][i]]]]
        assert child['conditions']['X_B'] == [original['conditions']['X_B'][i]]
        assert child['orchid_uncertainty']['values'][0]['budget'] == original['orchid_uncertainty']['values'][i]['budget']
        assert child['weight'] == round(500 / (500,1500)[i], 6)
        assert child['orchid_weighting']['parent_dataset'] == rel
        assert child['orchid_weighting']['measured_in'] == 'doi:shared'
        assert child['orchid']['observations'] == [{'id': f'ORC-{i}'}]
    before = copy.deepcopy(c['datasets']);assert E.split_unequal_uncertainties(c, pub) == []
    assert before == c['datasets']


def test_equal_sigmas_keep_group():
    c, pub, rel = fixture((.5,.5));before = copy.deepcopy(c['datasets'])
    assert E.split_unequal_uncertainties(c, pub) == [] and c['datasets'] == before


def test_solver_configurations_remain_aligned():
    c, pub, rel = fixture(solver=True); original = copy.deepcopy(c['datasets'][rel])
    E.split_unequal_uncertainties(c,pub)
    for i, child in enumerate(c['datasets'].values()):
        assert child['solver']['sublattice_configurations'] == [original['solver']['sublattice_configurations'][i]]
        assert child['solver']['sublattice_occupancies'] == [original['solver']['sublattice_occupancies'][i]]


def test_unaligned_arrays_fail_closed():
    c,pub,rel=fixture();c['datasets'][rel]['conditions']['X_B']=[.1]
    with pytest.raises(ValueError,match='unaligned'):E.split_unequal_uncertainties(c,pub)


def test_zpf_is_not_reweighted_from_temperature_uncertainty():
    c,pub,rel=fixture();c['datasets'][rel]['output']='ZPF';before=copy.deepcopy(c['datasets'])
    assert E.split_unequal_uncertainties(c,pub)==[] and c['datasets']==before
