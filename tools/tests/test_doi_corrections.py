#!/usr/bin/env python3
"""DOI correction table (tools/doi_corrections.json) and its scoped application.

Unit part (always): the table validates; lookups honour the context and scope (10.2320/jinstmet.73.180 and
10.1016/j.tca.2014.01.002 are dropped only for their one TDB header, 10.1557/jmr everywhere, handbook entries only for
their system); the literature exporter's helpers null a dropped DOI with doi_status / doi_correction and keep the
printed citation. Release part (when the orchid-data checkout is present): no dropped DOI is left in its context.
usage: python3 tools/tests/test_doi_corrections.py"""
import gzip, json, os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import doi_corrections as D  # noqa: E402
import export_literature_data as E  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


C = D.load()
check("table validates", D.validate(C) == [])
check("22 entries, 19 DOIs", len(C.entries) == 22 and len({e["doi"] for e in C.entries}) == 19)
check("version present", bool(C.version))
# scope: context_only TDB-header entries
check("tca dropped for NiTi-14Hu.tdb", (C.tdb("NiTi-14Hu.tdb", "10.1016/j.tca.2014.01.002") or {}).get("id") == "doicorr-06")
check("tca kept for its own paper (literature)", C.literature("10.1016/j.tca.2014.01.002") is None)
check("tca kept for another TDB", C.tdb("NiTi-13Pov.tdb", "10.1016/j.tca.2014.01.002") is None)
check("jinstmet dropped for IrPt-09Yam.tdb", C.tdb("databases/sgte-binary/IrPt-09Yam.tdb", "10.2320/JINSTMET.73.180") is not None)
check("jinstmet kept for MoB-09Yam.tdb", C.tdb("MoB-09Yam.tdb", "10.2320/jinstmet.73.180") is None)
check("jinstmet kept as a cited work", C.literature("10.2320/jinstmet.73.180") is None)
# scope: everywhere
for ctx, e in (("literature", C.literature("10.1557/JMR")), ("unrelated tdb", C.tdb("XyZ-99Abc.tdb", "10.1557/jmr")),
               ("handbook", C.handbook("AU-SN", "https://doi.org/10.1557/jmr")), ("anywhere", C.anywhere("10.1557/jmr"))):
    check("10.1557/jmr dropped everywhere (%s)" % ctx, e is not None and e["scope"] == "everywhere")
check("anywhere() ignores context_only entries", C.anywhere("10.1016/j.tca.2014.01.002") is None)
# handbook: system scoped, element order free
check("handbook IN-TE dropped", C.handbook(["Te", "In"], "10.1007/bf02418464") is not None)
check("handbook other system kept", C.handbook("AG-PU", "10.1007/bf02418464") is None)
# mark()
o = C.mark({"doi": "10.1037/013887", "citation": "printed"}, C.by_id["doicorr-08"])
check("mark nulls the DOI, keeps the citation", o == {"doi": None, "citation": "printed", "doi_status": "dropped_mismatch",
                                                     "doi_correction": "doicorr-08"})

# exporter helpers
store = E.Store("/nonexistent", "/nonexistent")
lst = E.cited_list([{"label": "12", "text": "J. E. Kittl, T.B. Massalski, J. Inst. Met. 93 (1964-1965) 182-188.",
                     "doi": "10.1001/jama.1964.03060380116096"},
                    {"label": "13", "text": "K. Yamada et al.", "doi": "10.2320/jinstmet.73.180"}], store)
check("cited source: dropped DOI null with status", lst[0].get("doi") is None and lst[0].get("doi_status") == "dropped_mismatch"
      and lst[0].get("doi_correction") == "doicorr-01" and lst[0]["reference"].startswith("J. E. Kittl"))
check("cited source: paper key = bibliography no-DOI id", lst[0]["paper"] == C.by_id["doicorr-01"]["replacement_id"])
check("cited source: correct DOI untouched", lst[1]["doi"] == "10.2320/jinstmet.73.180" and lst[1]["paper"] == "doi:10.2320/jinstmet.73.180")
pp = store.finalize_papers()
p0 = pp[lst[0]["paper"]]
check("papers.json: doi null, status, printed citation", p0.get("doi") is None and p0["doi_status"] == "dropped_mismatch"
      and p0["citation"].startswith("J. E. Kittl") and "doi_url" not in p0)


def model(tdb, doi):
    return E.model_link({"model_evidence": {"tdb": {"name": tdb}, "model_citation": {"doi": doi}}})


m = model("NiTi-14Hu.tdb", "10.1016/j.tca.2014.01.002")
check("model link: TDB-header DOI dropped", m.get("model_doi") is None and m.get("model_doi_status") == "dropped_mismatch"
      and m.get("model_doi_correction") == "doicorr-06")
check("model link: same DOI kept for another TDB", model("CrNiTi-14Hu.tdb", "10.1016/j.tca.2014.01.002").get("model_doi")
      == "10.1016/j.tca.2014.01.002")
E.SHA2DOI.clear()
E.SHA2DOI["bbbbbbbbbbbbbbbb"] = "10.1557/jmr.2008.0035"
check("paper key: journal-level DOI falls back to the PDF's own DOI",
      E.paper_key({"doi": "10.1557/JMR", "file_sha256": "b" * 64}) == "doi:10.1557/jmr.2008.0035")
check("paper key: journal-level DOI, unknown PDF -> sha256 key",
      E.paper_key({"doi": "10.1557/JMR", "file_sha256": "c" * 64}) == "sha256:" + "c" * 16)
a1 = {"id": "source:doi:10-1557-jmr", "doi": "10.1557/JMR", "file_sha256": "a" * 64}
a2 = {"id": "source:doi:10-1557-jmr", "doi": "10.1557/JMR", "file_sha256": "c" * 64}
store.add_source(a1)
store.add_source(a2)
check("colliding journal-level source ids: the file's own source wins",
      store.source_in_file("source:doi:10-1557-jmr", {"source:doi:10-1557-jmr": a2}) is a2)
check("journal-level DOI never maps a PDF", "a" * 16 not in E.SHA2DOI)

# release part
DATA = os.environ.get("ORCHID_DATA") or os.path.join(os.path.dirname(E.SITE), "orchid-data")
if os.path.isdir(os.path.join(DATA, "literature")):
    dropped_lit = {e["doi"] for e in C.entries if e["context"] == "literature_cited_source"}
    papers = json.load(open(os.path.join(DATA, "literature", "papers.json"), encoding="utf-8"))
    check("release: papers.json carries no literature-dropped DOI",
          not any(p.get("doi") in dropped_lit or k[4:] in dropped_lit for k, p in papers.items()))
    bad = 0
    for line in gzip.open(os.path.join(DATA, "literature", "orchid-literature.jsonl.gz"), "rt", encoding="utf-8"):
        r = json.loads(line)
        m = r.get("model") or {}
        if m.get("model_doi") and C.tdb(m.get("tdb"), m["model_doi"]):
            bad += 1
        if any(c.get("doi") in dropped_lit for c in r.get("cited_sources") or []):
            bad += 1
        if (r.get("paper") or "")[4:] in dropped_lit:
            bad += 1
    check("release: no dropped DOI in literature values (paper, cited sources, model links)", bad == 0)
    pd = json.load(open(os.path.join(DATA, "phase-diagrams", "sgte", "NiTi-14Hu.json"), encoding="utf-8"))["original_source"]
    check("release: phase-diagram NiTi-14Hu source DOI dropped, citation kept", pd.get("doi") is None and pd.get("citation")
          and pd.get("doi_correction") == "doicorr-06")
    ap = os.path.join(DATA, "assessed-phase-assemblages", "records.jsonl.gz")
    if os.path.exists(ap):
        n = 0
        for line in gzip.open(ap, "rt", encoding="utf-8"):
            r = json.loads(line)
            for c in r["citations"]:
                if c.get("doi") and C.handbook(r["system"], c["doi"]):
                    n += 1
        check("release: no handbook-dropped DOI in assessed assemblages", n == 0)
    bc = os.path.join(DATA, "bibliography", "doi_corrections.json")
    check("release: bibliography publishes the table", os.path.exists(bc) and open(bc, "rb").read() == open(D.TABLE_PATH, "rb").read())
else:
    print("SKIP release checks (no orchid-data checkout at %s)" % DATA)

print("%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
