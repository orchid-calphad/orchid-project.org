"""Describe stored boundary temperatures without changing the calculation domain."""
import math

def temperature_extent(record):
    values = [t for region in record.get('regions', []) for branch in region
              for t in branch.get('T', [])]
    if not values or not all(math.isfinite(t) and t > 0 for t in values):
        raise ValueError('Boundary temperatures must be nonempty, finite and positive')
    requested = record['T_range']
    if len(requested) != 2 or not all(math.isfinite(t) and t > 0 for t in requested) or requested[0] >= requested[1]:
        raise ValueError('Invalid requested temperature range')
    extent = [min(values), max(values)]
    return {'requested_K': list(requested), 'boundary_K': extent,
            'outside_requested_range': extent[0] < requested[0] or extent[1] > requested[1],
            'derivation': 'min/max of stored regions[*][*].T; no interpolation or clipping',
            'interpretation': 'Stored output extent, not a validated temperature domain or proof of equilibrium accuracy'}
