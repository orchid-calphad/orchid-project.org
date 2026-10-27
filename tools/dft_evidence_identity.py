"""Source-verified identity aliases; retain all literature values, count once in DFT fits."""
REVIEW = {
    'source_doi': '10.1016/j.jallcom.2009.09.036',
    'printed_page': 148, 'pdf_page': 3, 'table': 1,
    'basis': 'Visual inspection of Table 1: one This work formation enthalpy for each compound; duplicated extraction records point to the same printed result.',
    'pairs': [
        {'canonical_id': 'ORC-3UIVN2KBF5', 'alias_id': 'ORC-RN3WRMNFC6', 'compound': 'Ag2Dy', 'value_kJ_mol_atom': -32.579},
        {'canonical_id': 'ORC-KVFIPFQS3B', 'alias_id': 'ORC-P7PKLGDQ54', 'compound': 'Ag2Er', 'value_kJ_mol_atom': -34.119},
    ],
    'rule': 'Preserve both literature records and their numerical values; use the canonical identity once in DFT calibration and fitting exports.'
}

REVIEW["additional_reviews"] = [{'source_doi': '10.1016/j.tca.2013.05.024',
  'printed_page': 47,
  'pdf_page': 4,
  'table': 4,
  'basis': 'The 2022 assessment Table 2 (printed page 130, PDF page 5) repeats two first-principles '
           'formation enthalpies from Wang et al. 2013 Table 4 and cites that paper as reference 16; '
           'the reference text is on PDF page 11.',
  'reporting_doi': '10.1007/s11669-022-00947-8',
  'reporting_pdf_page': 5,
  'reporting_table': 2,
  'reporting_reference_label': '16',
  'source_hashes': {'10.1016/j.tca.2013.05.024': '46f3c7f049257145c80c983e9fc48df7f2c6f76ba0f6491c9eacf1337e28108f',
                    '10.1007/s11669-022-00947-8': 'fa9c11d530e12e0953b4900bb9f78dfef9400b6e9f2a1640bc7ee7a2119e6d5b'},
  'pairs': [{'canonical_id': 'ORC-6RCIZWHPCE',
             'alias_id': 'ORC-PDVG7NZMQP',
             'compound': 'Bi3Er5',
             'value_kJ_mol_atom': -66.035,
             'record_evidence': {'ORC-6RCIZWHPCE': {'paper': 'doi:10.1016/j.tca.2013.05.024',
                                                    'phase_label': 'Bi3Er5',
                                                    'temperature_K': 0.0,
                                                    'reference_state': {'BI': {'phase': 'RHOMBOHEDRAL_A7'},
                                                                        'ER': {'phase': 'HCP_A3'}}},
                                 'ORC-PDVG7NZMQP': {'paper': 'doi:10.1007/s11669-022-00947-8',
                                                    'phase_label': 'Er5Bi3',
                                                    'temperature_K': 0.0,
                                                    'reference_state': {'BI': {'phase': 'RHOMBOHEDRAL_A7'},
                                                                        'ER': {'phase': 'HCP_A3'}}}}},
            {'canonical_id': 'ORC-CWNSNQEM3C',
             'alias_id': 'ORC-6Y27WU6GQL',
             'compound': 'BiEr',
             'value_kJ_mol_atom': -82.069,
             'record_evidence': {'ORC-CWNSNQEM3C': {'paper': 'doi:10.1016/j.tca.2013.05.024',
                                                    'phase_label': 'BiEr',
                                                    'temperature_K': 0.0,
                                                    'reference_state': {'BI': {'phase': 'RHOMBOHEDRAL_A7'},
                                                                        'ER': {'phase': 'HCP_A3'}}},
                                 'ORC-6Y27WU6GQL': {'paper': 'doi:10.1007/s11669-022-00947-8',
                                                    'phase_label': 'α-ErBi',
                                                    'temperature_K': 0.0,
                                                    'reference_state': {'BI': {'phase': 'RHOMBOHEDRAL_A7'},
                                                                        'ER': {'phase': 'HCP_A3'}}}}}],
  'source_link_discrepancy': {'existing_resolved_doi': '10.1007/s11669-013-0231-4',
                              'verified_original_doi': '10.1016/j.tca.2013.05.024',
                              'basis': 'The existing DOI resolves to a La-B/La-Bi article, whereas the '
                                       'printed reference names Wang et al., Thermochimica Acta 566 '
                                       '(2013) 44-49. Identity is established from the printed '
                                       'reference and original PDF, not the erroneous resolved link.',
                              'status': 'corrected',
                              'correction': 'bier-reference-16-20261026'},
  'rule': 'Preserve both literature records; use the original 2013 calculation once in calibration and '
          'fitting.'}]

REVIEW["additional_reviews"].append({'source_doi': '10.1016/j.calphad.2013.07.013',
 'printed_page': 41,
 'pdf_page': 4,
 'table': 4,
 'reporting_doi': '10.1016/j.jallcom.2016.05.116',
 'reporting_pdf_page': 3,
 'reporting_table': 2,
 'reporting_reference_label': '43 (attribution in text); 46-49 (methods cited '
                              'in table)',
 'basis': 'Dong et al. 2013 Table 4 reports five ab initio Au-Ce formation '
          'enthalpies. Moustaine and Mahdouk 2016 Table 2 repeats them; the '
          'text above that table attributes the ab initio calculations to Dong '
          'et al. [43]. The numerical values and corresponding polymorphs '
          'agree.',
 'source_context_discrepancy': {'original_temperature_K': 0,
                                'later_report_temperature_K': 298.15,
                                'basis': 'The original Table 4 caption '
                                         'identifies ab initio results at 0 K. '
                                         'The later Table 2 prints 298.15 K '
                                         'for these repeated results. Both '
                                         'reported conditions and '
                                         'reference-state metadata remain in '
                                         'their records; the original '
                                         'calculation is canonical. Methods '
                                         'references 46-49 are preserved as '
                                         'printed, while reference 43 in the '
                                         'surrounding text establishes the '
                                         'calculation source.',
                                'status': 'documented; original source used '
                                          'for fitting'},
 'rule': 'Preserve both reports and their values and conditions; count the '
         'original calculation once in calibration and fitting.',
 'source_hashes': {'10.1016/j.calphad.2013.07.013': '4164a73884eaedc7eb52f6e282340557f2fa52cc2f5d1d8b01e8d3d580a0119a',
                   '10.1016/j.jallcom.2016.05.116': '0dff1db9f8641c1567f70fb1e397ce54a60ceb421f0f352cc3d213d08628867c'},
 'pairs': [{'canonical_id': 'ORC-G3N4GWDCFE.3a17a5',
            'alias_id': 'ORC-I5RYCRWVMX',
            'compound': 'AuCe2',
            'value_kJ_mol_atom': -50.99,
            'record_evidence': {'ORC-G3N4GWDCFE.3a17a5': {'paper': 'doi:10.1016/j.calphad.2013.07.013',
                                                          'phase_label': 'AuCe2',
                                                          'temperature_K': 0.0,
                                                          'reference_state': {'AU': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'},
                                                                              'CE': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'}}},
                                'ORC-I5RYCRWVMX': {'paper': 'doi:10.1016/j.jallcom.2016.05.116',
                                                   'phase_label': 'Ce2Au',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'AU': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'},
                                                                       'CE': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-JTTISQ65F3.427862',
            'alias_id': 'ORC-KI6DB4ZQSO',
            'compound': 'αAuCe',
            'value_kJ_mol_atom': -75.36,
            'record_evidence': {'ORC-JTTISQ65F3.427862': {'paper': 'doi:10.1016/j.calphad.2013.07.013',
                                                          'phase_label': 'αAuCe',
                                                          'temperature_K': 0.0,
                                                          'reference_state': {'AU': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'},
                                                                              'CE': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'}}},
                                'ORC-KI6DB4ZQSO': {'paper': 'doi:10.1016/j.jallcom.2016.05.116',
                                                   'phase_label': 'αCeAu',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'AU': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'},
                                                                       'CE': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-XFUVQWJIG7.eb4a61',
            'alias_id': 'ORC-PJWX7ZZH5C',
            'compound': 'Au2Ce',
            'value_kJ_mol_atom': -71.46,
            'record_evidence': {'ORC-XFUVQWJIG7.eb4a61': {'paper': 'doi:10.1016/j.calphad.2013.07.013',
                                                          'phase_label': 'Au2Ce',
                                                          'temperature_K': 0.0,
                                                          'reference_state': {'AU': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'},
                                                                              'CE': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'}}},
                                'ORC-PJWX7ZZH5C': {'paper': 'doi:10.1016/j.jallcom.2016.05.116',
                                                   'phase_label': 'CeAu2',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'AU': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'},
                                                                       'CE': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-HCS5Y5VHPZ.ad39d3',
            'alias_id': 'ORC-Y4KGJYHW7Y',
            'compound': 'Au6Ce',
            'value_kJ_mol_atom': -39.76,
            'record_evidence': {'ORC-HCS5Y5VHPZ.ad39d3': {'paper': 'doi:10.1016/j.calphad.2013.07.013',
                                                          'phase_label': 'Au6Ce',
                                                          'temperature_K': 0.0,
                                                          'reference_state': {'AU': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'},
                                                                              'CE': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'}}},
                                'ORC-Y4KGJYHW7Y': {'paper': 'doi:10.1016/j.jallcom.2016.05.116',
                                                   'phase_label': 'CeAu6',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'AU': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'},
                                                                       'CE': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-LRZ7JGSTYI.f7e67c',
            'alias_id': 'ORC-VDSEVXOHKO',
            'compound': 'βAuCe',
            'value_kJ_mol_atom': -74.71,
            'record_evidence': {'ORC-LRZ7JGSTYI.f7e67c': {'paper': 'doi:10.1016/j.calphad.2013.07.013',
                                                          'phase_label': 'βAuCe',
                                                          'temperature_K': 0.0,
                                                          'reference_state': {'AU': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'},
                                                                              'CE': {'phase': 'printed:Au '
                                                                                              '(Fcc), '
                                                                                              'Ce '
                                                                                              '(Fcc), '
                                                                                              'and '
                                                                                              'Sn '
                                                                                              '(Bct)'}}},
                                'ORC-VDSEVXOHKO': {'paper': 'doi:10.1016/j.jallcom.2016.05.116',
                                                   'phase_label': 'βCeAu',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'AU': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'},
                                                                       'CE': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'FCC_A1'}}}}}]})

REVIEW["additional_reviews"].append({'source_doi': '10.1016/j.intermet.2009.03.017',
 'printed_page': 881,
 'pdf_page': 4,
 'table': 5,
 'reporting_doi': '10.1016/j.calphad.2012.02.001',
 'reporting_pdf_page': 7,
 'reporting_table': 3,
 'reporting_reference_label': '55',
 'canonical_reporting_source': {'doi': '10.1016/j.calphad.2009.08.004',
                                'pdf_page': 3,
                                'printed_page': 675,
                                'table': 3,
                                'reference_label': '54'},
 'basis': 'Yuan et al. 2009 Table 3 and Schick et al. 2012 Table 3 both report '
          '-17.7 kJ/mol-atom for Mg2Si and cite Zhang et al. 2009. The '
          'original Zhang Table 5 confirms -17.70 kJ/mol-atom in the '
          'calculated column. The earlier reporting record is canonical among '
          'the two retained records; it is not a newly extracted '
          'original-source record.',
 'source_context_discrepancy': {'original_computation_temperature_K': 0,
                                'canonical_report_temperature_K': None,
                                'later_report_temperature_K': 298.15,
                                'basis': 'The original methods state that '
                                         'energies are calculated at 0 K. The '
                                         '2009 reporting record has no '
                                         'numerical temperature; the 2012 '
                                         'table labels its column Hf,298 '
                                         '(stored as 298.15 K under the '
                                         'existing documented convention). '
                                         'Preserve all reported temperatures '
                                         'and reference-state metadata, '
                                         'without adding a temperature to the '
                                         'earlier record.',
                                'status': 'documented; one underlying '
                                          'calculation counted'},
 'rule': 'Retain both literature reports, values and conditions. Count the '
         'earlier reporting record once in calibration and fitting, with the '
         'original calculation source explicitly identified.',
 'source_hashes': {'10.1016/j.intermet.2009.03.017': 'cbc333cfe7cd114f841317bb5cc315f61a0d33cb9ad6bba559f21cebd65e8cc6',
                   '10.1016/j.calphad.2009.08.004': '6ce3ec6f87678d37ab87b2f22c9f614a632cdbe5881a2ae52c853107c7bade06',
                   '10.1016/j.calphad.2012.02.001': 'f6d4439a1276ef71b8d53ec98750aa8d1c3902d168549f32daefd1c8ecffbe82'},
 'pairs': [{'canonical_id': 'ORC-AP4XHDIAMA',
            'alias_id': 'ORC-FYZFTH4VM6',
            'compound': 'Mg2Si',
            'value_kJ_mol_atom': -17.7,
            'record_evidence': {'ORC-AP4XHDIAMA': {'paper': 'doi:10.1016/j.calphad.2009.08.004',
                                                   'phase_label': 'Mg2Si',
                                                   'temperature_K': None,
                                                   'reference_state': {'MG': {'phase': 'HCP_A3'},
                                                                       'SI': {'phase': 'DIAMOND_A4'}}},
                                'ORC-FYZFTH4VM6': {'paper': 'doi:10.1016/j.calphad.2012.02.001',
                                                   'phase_label': 'Mg2Si',
                                                   'temperature_K': 298.15,
                                                   'reference_state': {'MG': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'HCP_A3'},
                                                                       'SI': {'fixed_state_variables': {'T': 298.15},
                                                                              'phase': 'DIAMOND_A4'}}}}}]})

REVIEW["additional_reviews"].append({'source_doi': '10.1016/j.jallcom.2009.05.083',
 'printed_page': 827,
 'pdf_page': 6,
 'table': 4,
 'reporting_doi': '10.1016/j.calphad.2011.09.007',
 'reporting_pdf_page': 2,
 'reporting_table': 2,
 'reporting_reference_label': '17',
 'canonical_reporting_source': {'doi': '10.1016/j.calphad.2009.07.001',
                                'pdf_page': 4,
                                'printed_page': 553,
                                'table': 2,
                                'reference_label': '29 (related '
                                                   'first-principles paper)'},
 'basis': 'Garay et al. CALPHAD 2009 Table 2 reports -46.9 kJ/mol-atom for '
          'SrSi in the first-principles GGA column. Li et al. 2011 Table 2 '
          'repeats -46.9 as GGA and cites Garay et al. J. Alloys Compd. 484 '
          '(2009) 822-831. That original Table 4 reports -46.89 in the GGA '
          'column, rounding to -46.9; the CALPHAD paper also cites that '
          'calculation paper as reference 29. The two retained reports refer '
          'to the same SrSi GGA result.',
 'source_rounding': {'original_source_value_kJ_mol_atom': -46.89,
                     'retained_reporting_value_kJ_mol_atom': -46.9,
                     'basis': 'Original Table 4 GGA column versus one-decimal '
                              'reports in the two CALPHAD tables. The retained '
                              'values are not replaced by the more precise '
                              'primary-source value.'},
 'source_link_discrepancy': {'existing_resolved_doi': '10.1016/j.calphad.2009.07.001',
                             'verified_original_doi': '10.1016/j.jallcom.2009.05.083',
                             'basis': 'Reference 17 in Li et al. names the J. '
                                      'Alloys Compd. calculation paper, but '
                                      'the current resolver links to the '
                                      'related CALPHAD assessment. Identity '
                                      'was verified against both source PDFs.',
                             'status': 'documented; broader '
                                       'reference-resolution correction '
                                       'remains separate'},
 'rule': 'Preserve both reporting records, values, uncertainties and '
         'conditions. Count the earlier available SrSi GGA report once in '
         'calibration and fitting. Keep LDA, PBE and other authors '
         'calculations distinct.',
 'source_hashes': {'10.1016/j.jallcom.2009.05.083': 'da47d3588297d938c1a092069456c39a4d9aad8ac56c6257773c731472114f0a',
                   '10.1016/j.calphad.2009.07.001': '146f5bb7dc5ec5f3b2b7ade53551be59be714306e16a4b42cad8678783dfb084',
                   '10.1016/j.calphad.2011.09.007': '0735078f7a34cd0438aca19b2566d871ba0d0e3e7c2a862c87f92c2fe5c8d48f'},
 'pairs': [{'canonical_id': 'ORC-5VW6EUBSRK.174c2f',
            'alias_id': 'ORC-ZFGCMVYTPX.bdc602',
            'compound': 'SrSi',
            'value_kJ_mol_atom': -46.9,
            'record_evidence': {'ORC-5VW6EUBSRK.174c2f': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'SrSi',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-ZFGCMVYTPX.bdc602': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'SiSr',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}}]})

REVIEW["additional_reviews"][3]["source_link_discrepancy"].update(status="corrected", correction="srsi-table2-20261026")
REVIEW["additional_reviews"].append({'source_doi': '10.1016/j.jallcom.2009.05.083',
 'printed_page': 827,
 'pdf_page': 6,
 'table': 4,
 'reporting_doi': '10.1016/j.calphad.2011.09.007',
 'reporting_pdf_page': 2,
 'reporting_table': 2,
 'reporting_reference_label': '17',
 'canonical_reporting_source': {'doi': '10.1016/j.calphad.2009.07.001',
                                'pdf_page': 4,
                                'printed_page': 553,
                                'table': 2,
                                'reference_label': '29 (related first-principles paper)'},
 'basis': 'Li et al. 2011 Table 2 repeats the same compound-specific GGA/LDA values reported in '
          'Garay et al. CALPHAD 2009 Table 2; the method and phase match. Preserve both reports '
          'and count the earlier available reporting record once.',
 'source_link_discrepancy': {'existing_resolved_doi': '10.1016/j.calphad.2009.07.001',
                             'verified_original_doi': '10.1016/j.jallcom.2009.05.083',
                             'basis': 'Reference 17 in Li et al. names the J. Alloys Compd. '
                                      'calculation paper, but the current resolver links to the '
                                      'related CALPHAD assessment. Identity was verified against '
                                      'both source PDFs.',
                             'status': 'corrected',
                             'correction': 'srsi-table2-20261026'},
 'rule': 'Preserve both reporting records, values, uncertainties and conditions. Count the earlier '
         'available SrSi GGA report once in calibration and fitting. Keep LDA, PBE and other '
         'authors calculations distinct.',
 'source_hashes': {'10.1016/j.jallcom.2009.05.083': 'da47d3588297d938c1a092069456c39a4d9aad8ac56c6257773c731472114f0a',
                   '10.1016/j.calphad.2009.07.001': '146f5bb7dc5ec5f3b2b7ade53551be59be714306e16a4b42cad8678783dfb084',
                   '10.1016/j.calphad.2011.09.007': '0735078f7a34cd0438aca19b2566d871ba0d0e3e7c2a862c87f92c2fe5c8d48f'},
 'pairs': [{'canonical_id': 'ORC-5BVWOE3YM2.88c8cb',
            'alias_id': 'ORC-XJH7MKZWNI.5dc3eb',
            'compound': 'Sr2Si',
            'value_kJ_mol_atom': -37.9,
            'record_evidence': {'ORC-5BVWOE3YM2.88c8cb': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'Sr2Si',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-XJH7MKZWNI.5dc3eb': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'SiSr2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-C4KDRUOFUL.8881a7',
            'alias_id': 'ORC-R7MCVZ3TXI.30dace',
            'compound': 'Sr2Si',
            'value_kJ_mol_atom': -41.3,
            'record_evidence': {'ORC-C4KDRUOFUL.8881a7': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'Sr2Si',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-R7MCVZ3TXI.30dace': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'SiSr2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-VELCVCLWDK.7808b8',
            'alias_id': 'ORC-4IQIRR4G6L.18917c',
            'compound': 'Sr5Si3',
            'value_kJ_mol_atom': -40.3,
            'record_evidence': {'ORC-VELCVCLWDK.7808b8': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'Sr5Si3',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-4IQIRR4G6L.18917c': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'Si3Sr5',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-YNLI4WOAVZ.3803d2',
            'alias_id': 'ORC-4NZW6V5CB7.2c7067',
            'compound': 'Sr5Si3',
            'value_kJ_mol_atom': -44.3,
            'record_evidence': {'ORC-YNLI4WOAVZ.3803d2': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'Sr5Si3',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-4NZW6V5CB7.2c7067': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'Si3Sr5',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-R5R6OJDLAW.254b80',
            'alias_id': 'ORC-VDSFOJ7ZWT.428991',
            'compound': 'SrSi',
            'value_kJ_mol_atom': -53.1,
            'record_evidence': {'ORC-R5R6OJDLAW.254b80': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'SrSi',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-VDSFOJ7ZWT.428991': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'SiSr',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-GOJ3G4HACO.c6ed69',
            'alias_id': 'ORC-UDVK5MDZST.136b26',
            'compound': 'β-SrSi2',
            'value_kJ_mol_atom': -35.9,
            'record_evidence': {'ORC-GOJ3G4HACO.c6ed69': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'β-SrSi2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-UDVK5MDZST.136b26': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'β-Si2Sr',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-CDZ7HLXDBU.ca7693',
            'alias_id': 'ORC-D5KD6WENLG.84dafb',
            'compound': 'β-SrSi2',
            'value_kJ_mol_atom': -44.7,
            'record_evidence': {'ORC-CDZ7HLXDBU.ca7693': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'β-SrSi2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-D5KD6WENLG.84dafb': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'β-Si2Sr',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-VFQ3MRKMJK.bca241',
            'alias_id': 'ORC-PPFJUUN6ZW.607caf',
            'compound': 'α-SrSi2',
            'value_kJ_mol_atom': -35.7,
            'record_evidence': {'ORC-VFQ3MRKMJK.bca241': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'α-SrSi2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-PPFJUUN6ZW.607caf': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'α-Si2Sr',
                                                          'temperature_K': None,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}},
           {'canonical_id': 'ORC-NIY3C4KGUT.f632d2',
            'alias_id': 'ORC-GM7U7E26UE.427288',
            'compound': 'α-SrSi2',
            'value_kJ_mol_atom': -42.0,
            'record_evidence': {'ORC-NIY3C4KGUT.f632d2': {'paper': 'doi:10.1016/j.calphad.2009.07.001',
                                                          'phase_label': 'α-SrSi2',
                                                          'temperature_K': 298.15,
                                                          'reference_state': {'SI': {'phase': 'DIAMOND_A4'},
                                                                              'SR': {'phase': 'FCC_A1'}}},
                                'ORC-GM7U7E26UE.427288': {'paper': 'doi:10.1016/j.calphad.2011.09.007',
                                                          'phase_label': 'α-Si2Sr',
                                                          'temperature_K': None,
                                                          'reference_state': {'SI': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'DIAMOND_A4'},
                                                                              'SR': {'fixed_state_variables': {'T': 298.15},
                                                                                     'phase': 'FCC_A1'}}}}}]})

def annotate(records):
    records = list(records)
    lookup = {r['id']: r for r in records}
    for review in [REVIEW] + REVIEW.get('additional_reviews', []):
        for pair in review['pairs']:
            ids = [pair['canonical_id'], pair['alias_id']]
            if not any(i in lookup for i in ids):
                continue
            if not all(i in lookup for i in ids):
                raise ValueError('Incomplete verified DFT identity pair: '+str(ids))
            for rid in ids:
                r = lookup[rid]
                expected = pair.get('record_evidence', {}).get(rid)
                paper = expected['paper'] if expected else 'doi:'+review['source_doi']
                label = expected['phase_label'] if expected else pair['compound']
                if (r.get('paper') != paper or r.get('unit') != 'kJ/mol-atom'
                    or r.get('value') != pair['value_kJ_mol_atom'] or r.get('property') != 'HM_FORM'
                    or r.get('evidence_class') != 'dft'
                    or label not in [p.get('label') for p in r.get('phases', [])]
                    or (expected and any(r.get(k) != expected[k] for k in ['temperature_K','reference_state']))):
                    raise ValueError('Verified DFT identity evidence changed: '+rid)
                r['evidence_identity'] = {'canonical_id':pair['canonical_id'], 'alias_ids':[pair['alias_id']],
                    'role':'canonical' if rid == pair['canonical_id'] else 'duplicate_representation',
                    'basis':review['basis'], 'source_doi':review['source_doi'],
                    'printed_page':review['printed_page'], 'table':review['table']}
                if 'reporting_doi' in review:
                    r['evidence_identity']['reporting_source'] = {
                        'doi':review['reporting_doi'], 'pdf_page':review['reporting_pdf_page'],
                        'table':review['reporting_table'], 'reference_label':review['reporting_reference_label']}
                    for key in ('source_link_discrepancy', 'source_context_discrepancy', 'canonical_reporting_source', 'source_rounding'):
                        if key in review:
                            r['evidence_identity'][key] = review[key]
    return records

def canonical(record):
    return record.get('evidence_identity', {}).get('canonical_id', record['id']) == record['id']
