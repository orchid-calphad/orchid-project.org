import pathlib,sys,types
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
import phase_label_check as p
import phase_mapper as pm
import build_assessed_assemblages as a

def model(*sublattices):
 return {'models':[{'sublattices':[{'constituents':list(x)} for x in sublattices]}]}

def test_terminal_must_reach_element_endmember():
 assert p.terminal_mapping_error('(αZr)','GAZR_BG',['GA','ZR'],model(['GA'],['ZR']))
 assert p.terminal_mapping_error('(αZr)','HCP_A3',['GA','ZR'],model(['GA','ZR'],['VA'])) is None
 assert p.terminal_mapping_error('αGaZr','GAZR_BG',['GA','ZR'],model(['GA'],['ZR'])) is None

def test_polymorph_and_molecular_species_are_distinguished():
 assert p.terminal_mapping_error('(αAl)','DIAMOND_A4',['AL','TI'],model(['AL','TI']))
 assert p.terminal_mapping_error('(αSb)','DIAMOND_A4',['SB','SN'],model(['SB','SN']))
 assert p.terminal_mapping_error('(δ-Pu,Th)','TETRAGONAL_A6',['PU','TH'],model(['PU','TH']))
 assert p.terminal_mapping_error('H2','H2GAS',['H','TI'],model(['H2'])) is None

def test_compound_cannot_use_assessed_solution_fallback():
 rid='phase:s-ti:hcp-a3'
 fake=types.SimpleNamespace(PM=pm,system=lambda key:({('test.tdb','HCP_A3'):rid},{},{}),
  grids_of=lambda key:[{'tdb':'test.tdb','sha256':'checked'}],
  tdb=lambda name,axis:({'HCP_A3':{'pearson':'hP2','constituents':[['TI','S'],['VA']]}},{'HCP_A3':(0,1)},'checked'))
 assert a.Registry.map_label(fake,'S-TI','Ti1+xS',{'pearson_symbol':'hP2'},('S','TI'),[])[0] is None
 assert a.Registry.map_label(fake,'S-TI','ε',{'pearson_symbol':'hP2'},('S','TI'),[])[0]==rid

def test_ternary_terminal_does_not_inherit_binary_hydride_mapping():
 en=model(['CA'],['H']);en['structure']={'pearson':'oP12'}
 rule=types.SimpleNamespace(by_tdbname={'CA-H':{'CAH2_C37':en}})
 r={'id':'test-ca','system':'CA-H-LI','elements':['CA','H','LI'],
    'phases':[{'label':'(αCa)','model_phase':'CAH2_C37'}]}
 p.apply([r],rule)
 assert 'model_phase' not in r['phases'][0]
 assert r['phases'][0]['model_phase_rejected']['reason']=='terminal_element_absent_from_model_endmembers'

def test_compound_with_fcc_pearson_is_not_an_element_solution():
 rid='phase:pd-zr:fcc-a1'
 fake=types.SimpleNamespace(PM=pm,system=lambda key:({('test.tdb','FCC_A1'):rid},{},{}),
  grids_of=lambda key:[{'tdb':'test.tdb','sha256':'checked'}],
  tdb=lambda name,axis:({'FCC_A1':{'pearson':'cF4','constituents':[['PD','ZR'],['VA']]}},{'FCC_A1':(0,1)},'checked'))
 assert a.Registry.map_label(fake,'PD-ZR','γPdZr',{'pearson_symbol':'cF4','strukturbericht':'A1'},('PD','ZR'),[])[0] is None
