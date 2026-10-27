import hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from check_espei_fit import check

def put(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d))
def fixture(tmp_path):
 b=tmp_path/'espei';b.mkdir();tdb=tmp_path/'AB.tdb'
 tdb.write_text('ELEMENT AL FCC_A1 26.98 0 0 !\nELEMENT NI FCC_A1 58.69 0 0 !\nTYPE_DEFINITION % SEQ * !\nPHASE LIQUID % 1 1 !\nCONSTITUENT LIQUID :AL,NI: !\nPARAMETER G(LIQUID,AL;0) 298 0; 6000 N !\nPARAMETER G(LIQUID,NI;0) 298 0; 6000 N !\n')
 sha=hashlib.sha256(tdb.read_bytes()).hexdigest();path='experimental/AL-NI/AB.hmr.json'
 d={'components':['AL','NI'],'phases':['LIQUID'],'output':'HMR','conditions':{'T':1000,'P':101325,'X_NI':[.5]},'values':[[[-1000]]],'weight':.5,'orchid':{'observations':[{'id':'ORC-one'}]},'orchid_uncertainty':{'default_sigma':{'value':500},'values':[{'sigma_effective_used':1000,'weight_used':.5}]}}
 put(b/path,d);put(b/'index.json',{'datasets':[{'path':path,'tdb':'AB.tdb'}]});put(tmp_path/'phases/systems/AL-NI.json',{'phases':[{'id':'liquid','models':[{'tdb':'AB.tdb','phase_name':'LIQUID','sha256':sha}]}]})
 m={'target':{'path':'AB.tdb','sha256':sha},'datasets':[{'path':path,'sha256':hashlib.sha256((b/path).read_bytes()).hexdigest()}]};sel=tmp_path/'selection.json';put(sel,m);return b,sel,m,d

def test_valid_selection(tmp_path):
 b,s,m,d=fixture(tmp_path);assert check(b,s)['ready']
def test_hash_change_blocks(tmp_path):
 b,s,m,d=fixture(tmp_path);m['target']['sha256']='bad';put(s,m);assert not check(b,s)['ready']
def test_repeated_input_blocks(tmp_path):
 b,s,m,d=fixture(tmp_path);m['datasets']*=2;put(s,m);assert any(e['code']=='repeated_dataset_path' for e in check(b,s)['errors'])
def test_override_must_describe_actual_fit_file(tmp_path):
 b,s,m,d=fixture(tmp_path);d['weight']=.25;put(tmp_path/'fit.json',d);m['datasets'][0]['fit_file']='fit.json';put(s,m)
 assert any(e['code']=='undocumented_weight_override' for e in check(b,s)['errors'])
 m['datasets'][0]['weight_override']={'value':.25,'reason':'Documented sensitivity experiment'};put(s,m);assert check(b,s)['ready']
def test_changed_measurement_cannot_be_hidden_by_override(tmp_path):
 b,s,m,d=fixture(tmp_path);d['values']=[[[-2000]]];put(tmp_path/'fit.json',d);m['datasets'][0]['fit_file']='fit.json';put(s,m)
 assert any(e['code']=='fit_file_changed_beyond_weight' for e in check(b,s)['errors'])
def test_missing_phase_blocks(tmp_path):
 b,s,m,d=fixture(tmp_path);d['phases']=['UNKNOWN'];p=b/m['datasets'][0]['path'];put(p,d);m['datasets'][0]['sha256']=hashlib.sha256(p.read_bytes()).hexdigest();put(s,m)
 assert any(e['code']=='target_phases_missing' for e in check(b,s)['errors'])

def add_second(b,s,m,d):
 import copy
 d=copy.deepcopy(d);p='experimental/AL-NI/AB.other.json';put(b/p,d)
 idx=json.loads((b/'index.json').read_text());idx['datasets'].append({'path':p,'tdb':'AB.tdb'});put(b/'index.json',idx)
 m['datasets'].append({'path':p,'sha256':hashlib.sha256((b/p).read_bytes()).hexdigest()});put(s,m);return p

def test_duplicate_evidence_across_files_blocks(tmp_path):
 b,s,m,d=fixture(tmp_path);add_second(b,s,m,d)
 assert any(e['code']=='repeated_evidence' for e in check(b,s)['errors'])
def test_derived_overlap_requires_specific_acknowledgement(tmp_path):
 b,s,m,d=fixture(tmp_path);d['orchid']={'points':[{'source_observation':'ORC-one'}]};d['values']=[[[-1100]]];other=add_second(b,s,m,d)
 r=check(b,s);assert not r['ready'];assert any(e['code']=='shared_derived_evidence' for e in r['review_required'])
 m['acknowledgements']=[{'code':'shared_derived_evidence','paths':[m['datasets'][0]['path'],other],'reason':'Both constraints intentionally used; correlated, not independent evidence.'}];put(s,m);assert check(b,s)['ready']
def test_solver_site_ratio_mismatch_blocks(tmp_path):
 b,s,m,d=fixture(tmp_path);d['solver']={'mode':'manual','sublattice_site_ratios':[2],'sublattice_configurations':[['AL']]};p=b/m['datasets'][0]['path'];put(p,d);m['datasets'][0]['sha256']=hashlib.sha256(p.read_bytes()).hexdigest();put(s,m)
 assert any(e['code']=='solver_sublattice_mismatch' for e in check(b,s)['errors'])
