"""Regression: Celsius ranges and uncertain values retain their widths after a page correction."""
import copy
import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import literature_gates as gates

class Report:
    def add(self, *args, **kwargs):
        pass

def test_range_correction_preserves_width_and_is_idempotent():
    record = {'id': 'range', 'property': 'TRANSITION_TEMPERATURE', 'value': 477.5,
              'temperature_K': 477.5, 'temperature_range_K': {'min': 400., 'max': 555.},
              'temperature_uncertainty': {'basis': 'printed_range', 'value': 77.5}, 'flags': {}}
    visual = {'page_fix:range': {'verdict': 'corrected', 'fields': {
        'temperature_K': 750.65, 'temperature_range_K': {'min': 673.15, 'max': 828.15}}}}
    gates.g0_page_fixes([record], visual, Report())
    assert record['value'] == record['temperature_K'] == 750.65
    assert record['temperature_range_K'] == {'min': 673.15, 'max': 828.15}
    assert record['temperature_range_K_stored'] == {'min': 400., 'max': 555.}
    assert record['temperature_uncertainty']['value'] == 77.5
    after = copy.deepcopy(record)
    gates.g0_page_fixes([record], visual, Report())
    assert record == after

def test_uncertain_point_and_unchecked_record():
    record = {'id': 'point', 'property': 'TRANSITION_TEMPERATURE', 'value': 186.7,
              'temperature_K': 186.7, 'temperature_uncertainty': {'basis': 'printed', 'value': 2.4}, 'flags': {}}
    untouched = copy.deepcopy(record); untouched['id'] = 'unchecked'
    original = copy.deepcopy(untouched)
    gates.g0_page_fixes([record, untouched], {'page_fix:point': {'verdict': 'corrected',
        'fields': {'temperature_K': 459.85}}}, Report())
    assert record['value'] == record['temperature_K'] == 459.85
    assert record['temperature_uncertainty']['value'] == 2.4
    assert untouched == original

def test_page_withdrawn_temperature_is_not_reinferred():
    import literature_r2_thermo as r2
    import literature_r3_thermo as r3
    record = {'id': 'various', 'property': 'HM_FORM', 'value': -58.3,
              'temperature_K': 298.15, 'temperature_printed': '298',
              'property_printed': 'enthalpy at 298 K', 'flags': {}}
    gates.g0_page_fixes([record], {'page_fix:various': {'verdict': 'corrected',
        'fields': {'temperature_K': None, 'temperature_basis': 'page_checked_various'}}}, Report())
    r2.r2_temperature_in_label([record], Report())
    r3.r3_temperature_printed([record], {}, Report())
    assert record['temperature_K'] is None
    assert record['value'] == -58.3
    assert record['temperature_K_stored'] == 298.15

def test_stability_bound_preserved_separately_from_measurement():
    record = {'id': 'bound', 'property': 'LATTICE_C', 'value': 6.419,
              'temperature_K': 595., 'flags': {}}
    bound = {'op': '<', 'value_K': 868.15, 'printed': '<595 °C'}
    gates.g0_page_fixes([record], {'page_fix:bound': {'verdict': 'corrected',
        'fields': {'temperature_K': None, 'temperature_bound': bound}}}, Report())
    assert record['temperature_K'] is None
    assert record['temperature_bound'] == bound
    assert record['value'] == 6.419

def test_celsius_offset_does_not_invent_reading_precision():
    import uncertainty_budget as ub
    before = {'property': 'TRANSITION_TEMPERATURE', 'value': 186.7}
    after = {'property': 'TRANSITION_TEMPERATURE', 'value': 459.85,
             'corrections': [{'field': 'temperature_K', 'before': 186.7, 'after': 459.85,
                              'rule': 'printed_Celsius_to_Kelvin'}]}
    assert ub.reading_component(after, 'value') == ub.reading_component(before, 'value')
    assert ub.reading_component(after, 'value')['value'] == 0.05

def test_enthalpy_reclassification_preserves_numeric_value():
    record = {'id': 'enthalpy', 'property': 'INVARIANT_TEMPERATURE', 'value': 19000.,
              'temperature_K': 19000., 'flags': {}}
    gates.g0_page_fixes([record], {'page_fix:enthalpy': {'verdict': 'corrected',
        'fields': {'value': 19000., 'temperature_K': None,
                   'property': 'TRANSFORMATION_ENTHALPY', 'unit': 'J/mol'}}}, Report())
    assert record['value'] == 19000.
    assert record['temperature_K'] is None
    assert record['property'] == 'TRANSFORMATION_ENTHALPY'
