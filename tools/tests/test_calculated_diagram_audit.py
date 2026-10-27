import hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit_calculated_diagrams import structural,phase_name

def fixture(tmp_path):
    tdb=tmp_path/'AB.tdb';tdb.write_text('source fixture')
    d={'tdb':tdb.name,'sha256':hashlib.sha256(tdb.read_bytes()).hexdigest(),'phases':['LIQUID'],'T_range':[300,1200],
       'regions':[[{'phase':'LIQUID','x':[.1,.2,.3],'T':[500,600,700]},
                   {'phase':'LIQUID#2','x':[.9,.8,.7],'T':[500,600,700]}]],
       'invariants':[{'T':500,'points':[{'phase':'LIQUID','x':.1},{'phase':'LIQUID#2','x':.9}]}]}
    p=tmp_path/'diagram.json';p.write_text(json.dumps(d));return p,tdb,d

def test_composition_sets_and_invariant_avoidance(tmp_path):
    p,tdb,d=fixture(tmp_path);r=structural(p,tmp_path)
    assert not r['issues']
    assert r['probes'][0]['T_K']==600
    assert r['probes'][0]['expected']==['LIQUID','LIQUID']
    assert phase_name('PHASE#12')=='PHASE'
    assert phase_name('PHASE_NAME')=='PHASE_NAME'

def test_source_mismatch_and_invalid_coordinates(tmp_path):
    p,tdb,d=fixture(tmp_path);tdb.write_text('changed')
    d['regions'][0][0]['x'][0]=2;p.write_text(json.dumps(d))
    r=structural(p,tmp_path)
    assert 'source_hash_mismatch' in r['issues']
    assert 'invalid_boundary_coordinate' in r['issues']

def test_grid_run_lengths_are_checked(tmp_path):
    p,tdb,d=fixture(tmp_path)
    d['grid']={'x_points':3,'x_step':.5,'T0':300,'dT':10,'labels':['LIQUID'],'rows':[[[0,2]]]}
    p.write_text(json.dumps(d));assert 'grid_row_length_mismatch' in structural(p,tmp_path)['issues']

def test_empty_map_is_review_flag_not_proof_of_bad_thermodynamics(tmp_path):
    p,tdb,d=fixture(tmp_path);d['regions']=[];p.write_text(json.dumps(d))
    r=structural(p,tmp_path)
    assert 'empty_boundaries_requires_review' in r['warnings']
    assert 'no_checkable_probe' in r['warnings']


def test_documented_extent_does_not_hide_bad_coordinates(tmp_path):
    from diagram_temperature_extent import temperature_extent
    p,tdb,d=fixture(tmp_path)
    d['regions'][0][0]['T'][0]=299
    d['temperature_extent']=temperature_extent(d)
    p.write_text(json.dumps(d));r=structural(p,tmp_path)
    assert 'invalid_boundary_coordinate' not in r['issues']
    assert 'documented_boundary_extent_outside_requested_range' in r['warnings']
    d['temperature_extent']['boundary_K'][0]=298
    p.write_text(json.dumps(d));r=structural(p,tmp_path)
    assert 'temperature_extent_metadata_mismatch' in r['issues']
    assert 'invalid_boundary_coordinate' in r['issues']
    d['temperature_extent']=temperature_extent(d)
    d['regions'][0][0]['x'][0]=2
    p.write_text(json.dumps(d))
    assert 'invalid_boundary_coordinate' in structural(p,tmp_path)['issues']
