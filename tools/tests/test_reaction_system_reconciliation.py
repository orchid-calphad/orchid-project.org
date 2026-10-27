"""Preserve classification and fit safeguards after restoring reaction partners."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import literature_gates as gates


class ReactionSystemReconciliationTest(unittest.TestCase):
    def record(self):
        return {
            'id': 'reaction-system-regression', 'property': 'INVARIANT_TEMPERATURE',
            'system': 'CE', 'source_system': 'CA-CE', 'elements': ['CE'],
            'reaction': {'printed': 'Liquid 1 -> bcc(Ce) + bcc(Ca)'},
            'phases': [{'label': p} for p in ['Liquid 1', 'bcc(Ce)', 'bcc(Ca)']],
            'flags': {'espei_ready': False, 'espei_certified': False},
            'gate_failures': [{'code': 'phases_added_from_printed_reaction'}],
        }

    def test_restores_binary_class_without_promoting_fit(self):
        record = self.record()
        gates.reconcile_reaction_systems([record], gates.Report())
        self.assertEqual(record['elements'], ['CA', 'CE'])
        self.assertEqual(record['reaction']['type'], 'eutectic')
        self.assertEqual(record['reaction']['type_basis'], 'derived_from_printed_reaction')
        self.assertFalse(record['flags']['espei_ready'])
        self.assertFalse(record['flags']['espei_certified'])
        captured = copy.deepcopy(record)
        gates.reconcile_reaction_systems([record], gates.Report())
        self.assertEqual(record, captured)

    def test_preserves_recorded_class(self):
        record = self.record()
        record['reaction']['type'] = 'peritectic'
        gates.reconcile_reaction_systems([record], gates.Report())
        self.assertEqual(record['reaction']['type'], 'peritectic')
        self.assertFalse(record['flags']['espei_ready'])

    def test_does_not_expand_beyond_source_system(self):
        record = self.record()
        record['source_system'] = 'CE'
        captured = copy.deepcopy(record)
        gates.reconcile_reaction_systems([record], gates.Report())
        self.assertEqual(record, captured)


if __name__ == '__main__':
    unittest.main()
