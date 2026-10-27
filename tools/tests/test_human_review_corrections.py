"""Regression checks for guarded, atomic human-review correction batches."""
import copy
import json
from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
from human_review_corrections import apply


class ReviewCorrectionsTest(unittest.TestCase):
    def setUp(self):
        self.patches = json.loads((TOOLS/'human_review_corrections_20261026.json').read_text())
        self.records = [dict(copy.deepcopy(p['before']), id=p['record_id'],
                             flags={'espei_ready': True, 'espei_certified': True,
                                    'plausibility': {'result': 'pass'}}) for p in self.patches]

    def test_batch_and_repeat(self):
        result = apply(self.records, self.patches)
        self.assertEqual((len(result['applied']), result['fields_changed']), (19, 31))
        for record, patch in zip(self.records, self.patches):
            for key, value in patch['after'].items():
                self.assertEqual(record[key], value)
            self.assertFalse(record['flags']['espei_ready'])
            self.assertFalse(record['flags']['espei_certified'])
            self.assertNotIn('plausibility', record['flags'])
            self.assertEqual(record['plausibility_before_correction'], {'result': 'pass'})
        before = copy.deepcopy(self.records)
        self.assertEqual(apply(self.records, self.patches)['applied'], [])
        self.assertEqual(before, self.records)

    def test_conflict_aborts_every_record(self):
        self.records[-1]['paper'] = 'doi:conflicting-source'
        before = copy.deepcopy(self.records)
        with self.assertRaisesRegex(ValueError, 'changed guarded fields'):
            apply(self.records, self.patches)
        self.assertEqual(before, self.records)

    def test_missing_aborts_every_record(self):
        self.records.pop()
        before = copy.deepcopy(self.records)
        with self.assertRaisesRegex(ValueError, 'missing record'):
            apply(self.records, self.patches)
        self.assertEqual(before, self.records)

    def test_unrelated_record_preserved(self):
        extra = {'id': 'unrelated', 'value': 42, 'flags': {'espei_ready': True}}
        self.records.append(copy.deepcopy(extra))
        apply(self.records, self.patches)
        self.assertEqual(self.records[-1], extra)


if __name__ == '__main__':
    unittest.main()
