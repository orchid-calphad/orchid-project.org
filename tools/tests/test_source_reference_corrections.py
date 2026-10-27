import copy
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import source_reference_corrections as C


def records():
 rs=copy.deepcopy(C.REVIEW['records'])
 for r in rs:r['paper']='doi:'+C.REVIEW['reporting_doi']
 return rs


def test_all_reviewed_sources_corrected_without_altering_values():
 rs=records();before=copy.deepcopy(rs);C.annotate(rs)
 for a,b in zip(before,rs):
  assert {k:v for k,v in a.items() if k not in ('origin','cited_sources')}=={k:v for k,v in b.items() if k not in ('origin','cited_sources')}
  assert b['origin']['measured_in']=='doi:'+C.NEW
  assert b['origin']['measured_in_original']=='doi:'+C.OLD
  assert b['cited_sources'][0]['doi']==C.NEW
  assert b['cited_sources'][0]['reference']==a['cited_sources'][0]['reference']
 snap=copy.deepcopy(rs);C.annotate(rs);assert rs==snap


@pytest.mark.parametrize('field,value',[('value',1.234),('paper','doi:unrelated')])
def test_changed_evidence_rejected(field,value):
 rs=records();rs[0][field]=value
 with pytest.raises(ValueError):C.annotate(rs)


def test_other_papers_keep_valid_doi():
 rs=records();rs[0]['id']='unreviewed';rs=rs[:1];before=copy.deepcopy(rs)
 C.annotate(rs);assert rs==before


def test_different_printed_reference_rejected():
 rs=records();rs[0]['cited_sources'][0]['reference']='unrelated'
 with pytest.raises(ValueError):C.annotate(rs)


def test_paper_memberships_keep_the_unrelated_original_work():
 rs=records();C.annotate(rs)
 for r in rs:r['system']='BI-ER'
 rs.append({'id':'La-paper-value','paper':'doi:'+C.OLD,'system':'B-LA'})
 papers={'doi:'+C.OLD:{'roles':['reporting_paper','cited_source'],'source_kinds':['journal_article','cited_work']},'doi:'+C.NEW:{'roles':['reporting_paper'],'source_kinds':['journal_article']}}
 C.reconcile_papers(rs,papers)
 assert papers['doi:'+C.OLD]['systems']==['B-LA']
 assert papers['doi:'+C.OLD]['roles']==['reporting_paper']
 assert papers['doi:'+C.NEW]['roles']==['cited_source','reporting_paper']
