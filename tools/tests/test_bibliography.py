#!/usr/bin/env python3
"""Synthetic tests of tools/build_bibliography.py helpers: DOI normalisation, ids/keys, printed first-author parsing,
the Crossref parser (no abstract kept, licences, ORCIDs, reference DOIs), the open-access rule (Crossref licence
field only) and the DOI consistency check. No network, no data files.
usage: ~/orchid-py313/bin/python tools/tests/test_bibliography.py"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build_bibliography as B  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


check("norm_doi url", B.norm_doi("https://doi.org/10.1016/J.CALPHAD.2005.02.002") == "10.1016/j.calphad.2005.02.002")
check("norm_doi prefix", B.norm_doi("doi: 10.1007/BF00810049") == "10.1007/bf00810049")
check("norm_doi junk", B.norm_doi("not a doi") is None and B.norm_doi(None) is None)
check("key stable", B.key_of("doi:10.1016/j.calphad.2005.02.002") == "pd5e729a388c2c588")
a = B.nodoi_id("Kittl", 1964, "J. E. Kittl, T.B. Massalski, J. Inst. Met. 93")
check("nodoi id stable under case/punctuation", a == B.nodoi_id("KITTL", 1964, "j e kittl t b massalski j inst met 93"))
check("nodoi id prefix", a.startswith("nodoi:") and len(a) == 22)
for c, want in [("J. E. Kittl, T.B. Massalski, J. Inst. Met.", "kittl"), ("J. Alloys Compd. 385, 133", None),
                ("J. Thermal Anal. Calorim., 125", None), ("Z. Metallkd., 83, 21", None),
                ("A. Fernandez Guillermet, Metall. Trans.", "fernandez guillermet"),
                ("B. Joensson and J. Agren", "joensson"), ("Kittl, J.E., foo", "kittl")]:
    check(f"first author {c!r}", B.first_author_surname(c) == want)
rec = {"http_status": 200, "retrieved_at": "2026-10-23T12:00:00Z", "body": {"message": {
    "DOI": "10.1/X", "title": ["T"], "author": [{"given": "A.", "family": "Smith", "ORCID": "https://orcid.org/0000-0001-2345-6789",
                                                  "authenticated-orcid": True}],
    "issued": {"date-parts": [[1999, 5]]}, "type": "journal-article", "is-referenced-by-count": 7,
    "license": [{"URL": "https://www.elsevier.com/tdm/userlicense/1.0/", "content-version": "tdm",
                 "start": {"date-parts": [[1999, 5, 1]]}, "delay-in-days": 0},
                {"URL": "http://creativecommons.org/licenses/by/4.0/", "content-version": "vor",
                 "start": {"date-parts": [[2015, 1, 1]]}, "delay-in-days": 0}],
    "reference": [{"key": "r1", "DOI": "10.2000/Y", "doi-asserted-by": "publisher"}, {"key": "r2", "unstructured": "x"}],
    "abstract": "<jats:p>never kept</jats:p>"}}}
cr = B.parse_crossref(rec)
check("crossref fields", cr["year"] == 1999 and cr["authors"][0]["orcid"] == "0000-0001-2345-6789" and cr["is_referenced_by_count"] == 7)
check("crossref no abstract", "abstract" not in cr and "never kept" not in str(cr))
check("crossref references", cr["references_with_doi"] == [{"doi": "10.2000/y", "asserted_by": "publisher"}] and cr["references_deposited"] == 2)
oa = B.oa_license(cr["license"])
check("oa licence from crossref only", oa["url"] == "http://creativecommons.org/licenses/by/4.0/" and oa["content_version"] == "vor")
check("no oa without cc", B.oa_license([{"url": "https://www.elsevier.com/tdm/userlicense/1.0/"}]) is None)
check("404 -> None", B.parse_crossref({"http_status": 404, "body": None}) is None)
check("doi check consistent", B.doi_check({"first_author": "smith", "year": 1999}, cr)["status"] == "consistent")
check("doi check fuzzy (OCR)", B.doi_check({"first_author": "smlth", "year": 1999}, cr)["status"] == "consistent")
bad = B.doi_check({"first_author": "jones", "year": 1960}, cr)
check("doi check mismatch", bad["status"] == "mismatch_suspect" and set(bad["issues"]) == {"first_author_not_in_crossref_authors", "year_differs"})
check("system id", B.norm_system(["Ni", "Al"]) == "AL-NI" and B.norm_system("ti-ni") == "NI-TI")
print(f"{len(fails)} failures")
sys.exit(1 if fails else 0)
