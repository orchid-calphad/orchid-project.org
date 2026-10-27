import copy
from test_temperature_page_corrections import gates, Report

def test_composition_and_phase_identity_survive_page_correction():
 r={'id':'endmember','system':'NI-TI','elements':['NI','TI'],'property':'HM_FORM','value':7.08,'overall_composition':{'TI':.25},'phases':[{'label':'Ni3Ti'}],'flags':{}}
 old=copy.deepcopy(r)
 fields={'system':'NI','elements':['NI'],'overall_composition':{'NI':1.},'property':'END_MEMBER_FORMATION_ENERGY','phases':[{'label':'(Ni)3(Ni)1','parent_phase':'Ni3Ti'}]}
 visual={'page_fix:endmember':{'verdict':'corrected','fields':fields}}
 gates.g0_page_fixes([r],visual,Report())
 assert all(r[k]==v for k,v in fields.items())
 assert r['overall_composition_stored']==old['overall_composition']
 assert r['phases_stored']==old['phases']
 before=copy.deepcopy(r);gates.g0_page_fixes([r],visual,Report());assert r==before

def test_bound_withdraws_exact_value_but_retains_limit():
 r={'id':'bound','property':'INVARIANT_TEMPERATURE','temperature_K':393.,'value':393.,'flags':{}}
 fields={'value':None,'temperature_K':None,'temperature_bound':{'op':'<','value_K':393.,'printed':'<393'}}
 gates.g0_page_fixes([r],{'page_fix:bound':{'verdict':'corrected','fields':fields}},Report())
 assert r['value'] is None and r['temperature_K'] is None
 assert r['value_stored']==393. and r['temperature_bound']['value_K']==393.

def test_late_reaction_partners_expand_narrowed_system():
 r={'id':'agco','system':'AG','elements':['AG'],'source_system':'AG-CO','flags':{},
    'phases':[{'label':'fcc(Ag)'},{'label':'fcc(Co)'},{'label':'hcp(Co)'}],
    'gate_failures':[{'code':'phases_added_from_printed_reaction'}]}
 gates.reconcile_reaction_systems([r],Report())
 assert r['system']=='AG-CO' and r['elements']==['AG','CO']
 single={'id':'bismuth','system':'BI','elements':['BI'],'source_system':'AU-BI','flags':{},
    'phases':[{'label':'Bi(liq)'},{'label':'Bi(rho)'}],
    'gate_failures':[{'code':'phases_added_from_printed_reaction'}]}
 gates.reconcile_reaction_systems([single],Report())
 assert single['system']=='BI'
