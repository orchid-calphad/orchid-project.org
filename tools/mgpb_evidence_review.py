"""Source-verified Mg2Pb evidence corrections; numerical literature fields are immutable."""
import copy

DOI = '10.1016/j.jallcom.2020.158396'
GATE = 'experimental_conditions_require_primary_source_review'
REVIEW = {'source_doi': '10.1016/j.jallcom.2020.158396',
 'source_pdf_sha256': 'af6653177295b67e24a9852206942313bd0219e89430bd6396fb4444619473c5',
 'pages': [3, 4],
 'table': 3,
 'decision': 'Five rows described as experimental on page 3; both unit columns of each row represent '
             'one observation. Preserve all values; exclude from fitting and calibration until original '
             'conditions are verified.',
 'canonical_selection': 'Use the J/mol atom column as the identity anchor because it has finer printed '
                        'resolution; do not overwrite the eV/atom representation.',
 'discrepancies': [{'row': 7,
                    'field': 'value',
                    'table_J_mol_atom': -5855,
                    'text_J_mol_atom': -5866,
                    'action': 'retain table value; flag disagreement without substitution'},
                   {'row': 10,
                    'field': 'reference',
                    'table_reference': '24',
                    'text_reference': '25',
                    'action': 'retain printed citation; flag disagreement without substitution'}],
 'pairs': [{'row': 7,
            'canonical_id': 'ORC-BRBEV2CMMI.7deb20',
            'alias_id': 'ORC-BRBEV2CMMI.f3b1dd',
            'records': {'ORC-BRBEV2CMMI.7deb20': {'value': -5.855,
                                                  'value_printed': -5855.0,
                                                  'unit_printed': 'J/mol atom'},
                        'ORC-BRBEV2CMMI.f3b1dd': {'value': -5.847011126,
                                                  'value_printed': -0.0606,
                                                  'unit_printed': 'eV/atom'}}},
           {'row': 8,
            'canonical_id': 'ORC-VASZEU4BTD.47e16d',
            'alias_id': 'ORC-VASZEU4BTD.11e7a0',
            'records': {'ORC-VASZEU4BTD.11e7a0': {'value': -5.499663931,
                                                  'value_printed': -0.057,
                                                  'unit_printed': 'eV/atom'},
                        'ORC-VASZEU4BTD.47e16d': {'value': -5.5,
                                                  'value_printed': -5500.0,
                                                  'unit_printed': 'J/mol atom'}}},
           {'row': 9,
            'canonical_id': 'ORC-6NMZMDWEIZ.6c25ee',
            'alias_id': 'ORC-6NMZMDWEIZ.e16ce5',
            'records': {'ORC-6NMZMDWEIZ.6c25ee': {'value': -5.855,
                                                  'value_printed': -5855.0,
                                                  'unit_printed': 'J/mol atom'},
                        'ORC-6NMZMDWEIZ.e16ce5': {'value': -5.847011126,
                                                  'value_printed': -0.0606,
                                                  'unit_printed': 'eV/atom'}}},
           {'row': 10,
            'canonical_id': 'ORC-UN5ZR3DH4R.a5aded',
            'alias_id': 'ORC-UN5ZR3DH4R.9ee812',
            'records': {'ORC-UN5ZR3DH4R.9ee812': {'value': -5.52860953,
                                                  'value_printed': -0.0573,
                                                  'unit_printed': 'eV/atom'},
                        'ORC-UN5ZR3DH4R.a5aded': {'value': -5.533,
                                                  'value_printed': -5533.0,
                                                  'unit_printed': 'J/mol atom'}}},
           {'row': 11,
            'canonical_id': 'ORC-XSC6EPABZQ.4d8b85',
            'alias_id': 'ORC-XSC6EPABZQ.56a13f',
            'records': {'ORC-XSC6EPABZQ.4d8b85': {'value': -5.9,
                                                  'value_printed': -5900.0,
                                                  'unit_printed': 'J/mol atom'},
                        'ORC-XSC6EPABZQ.56a13f': {'value': -5.895253793,
                                                  'value_printed': -0.0611,
                                                  'unit_printed': 'eV/atom'}}}]}


def annotate(records):
    records = list(records)
    lookup = {r['id']: r for r in records}
    for pair in REVIEW['pairs']:
        ids = [pair['canonical_id'], pair['alias_id']]
        if not any(rid in lookup for rid in ids):
            continue
        if not all(rid in lookup for rid in ids):
            raise ValueError('Incomplete Mg2Pb unit pair: '+str(ids))
        for rid in ids:
            r = lookup[rid]
            if (r.get('paper') != 'doi:'+DOI or r.get('property') != 'HM_FORM'
                or r.get('unit') != 'kJ/mol-atom'
                or any(r.get(k) != v for k,v in pair['records'][rid].items())
                or 'Mg2Pb' not in [p.get('label') for p in r.get('phases', [])]):
                raise ValueError('Mg2Pb source evidence changed: '+rid)
            identity = {'canonical_id':pair['canonical_id'], 'alias_ids':[pair['alias_id']],
                        'role':'canonical' if rid == pair['canonical_id'] else 'alternate_unit_representation',
                        'basis':'Two printed unit columns of the same Table 3 row; not independent evidence',
                        'source_doi':DOI, 'pdf_page':4, 'table':3, 'row':pair['row']}
            if r.get('source_evidence_review', {}).get('review_id') == 'mgpb-table3-20261026':
                if r.get('evidence_class') != 'experimental' or r.get('temperature_K') is not None or r.get('evidence_identity') != identity:
                    raise ValueError('Mg2Pb correction changed: '+rid)
                continue
            if r.get('evidence_class') != 'dft' or r.get('temperature_K') != 0:
                raise ValueError('Unexpected Mg2Pb pre-correction state: '+rid)
            before = {k:copy.deepcopy(r.get(k)) for k in ['evidence_class','origin','method','method_basis','temperature_K','temperature_basis','flags']}
            r['evidence_class'] = 'experimental'
            r['origin'] = {'class':'experimental','basis':'page_3_text_explicitly_describes_experimental_results','own_work_of_paper':False}
            r['method'] = 'experimental, method not specified in reporting paper'
            r['method_basis'] = 'page_3_results_text; original experimental method not inferred'
            r.setdefault('temperature_K_stored', r['temperature_K'])
            r.setdefault('temperature_basis_stored', r.get('temperature_basis'))
            r['temperature_K'] = None
            r['temperature_basis'] = 'measurement_temperature_requires_primary_source_review'
            r['evidence_identity'] = identity
            r['source_evidence_review'] = {'review_id':'mgpb-table3-20261026','source_doi':DOI,
                'pdf_pages':[3,4],'table':3,'row':pair['row'],'before':before,
                'fitting_status':'excluded_pending_primary_source_conditions',
                'temperature_note':'The stored 0 K comes from the mixed table heading; it does not establish experimental measurement temperature. The printed heading is retained.',
                'discrepancies':[d for d in REVIEW['discrepancies'] if d['row']==pair['row']]}
            flags = r.setdefault('flags', {})
            flags.update(espei_ready=False, espei_certified=False, state_fully_specified=False,
                         review_state='source_checked_conditions_pending')
            flags['specification_gaps'] = ['state:experimental_temperature','state:primary_source_conditions']
            r.setdefault('gate_failures', []).append({'code':GATE,'effect':'not_espei_ready',
                'reason':'Experimental rows in a mixed DFT table; original measurement conditions unresolved. Both unit representations remain excluded from calibration and fitting.'})
    return records
