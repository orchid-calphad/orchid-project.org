#!/usr/bin/env python3
"""Fixture test: _export.reader_agreement of the record store is published per value as `reader_agreement` +
`reader_agreement_label` (export_literature_data.load_tables), the CSV carries both columns, and the second Claude reader
id 'claude-opus-5.5#2' stays distinguishable in `readers`. No real corpus is read; nothing is written outside a temp dir.
usage: python3 tools/tests/test_reader_agreement.py"""
import json, os, sys, tempfile
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import export_literature_data as E  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


TMP = tempfile.mkdtemp(prefix="orchid-ra-")
os.makedirs(f"{TMP}/db-v2.1/general/records")
src = {"record_type": "source", "id": "source:doi:10.1000/x", "kind": "journal_article", "doi": "10.1000/x"}
art = {"record_type": "artifact", "id": "artifact:abc:p3:table-1", "source": "source:doi:10.1000/x", "label": "1",
       "caption": "Table 1. Enthalpies of formation", "page": 3}


def obs(row, cls, readings):
    o = {"record_type": "observation", "id": f"obs:artifact:abc:p3:table-1:row-{row}", "system": ["Al", "Ni"],
         "observation_type": "equilibrium_property", "quantity": "HM_FORM", "value": -41700.0 - row, "unit": "J/mol",
         "origin": "experimental", "readings": readings,
         "quality": {"status": {"extraction_accepted": cls == "multi_reader"}}, "_export": {}}
    if cls:
        o["_export"].update({"reader_agreement": cls, "reader_agreement_detail": {"row_state": "x", "agreement": "x", "visual_readers": []}})
    return o


r0 = obs(0, "multi_reader", ["reading:a:row-0:gemini-3-8-flash", "reading:a:row-0:claude-opus-5-5-2"])
r1 = obs(1, "single_reader", ["reading:a:row-1:gemini-3-8-flash"])
r2 = obs(2, "single_reader_text_layer", ["reading:a:row-2:claude-opus-5-5", "reading:a:row-2:text-layer"])
r3 = obs(3, None, ["reading:a:row-3:gemini-3-8-flash", "reading:a:row-3:claude-opus-5-5"])
with open(f"{TMP}/db-v2.1/general/records/fixture.jsonl", "w") as s:
    for r in (src, art, r0, r1, r2, r3):
        s.write(json.dumps(r) + "\n")

store = E.Store(TMP, TMP)
records = []
E.load_tables(store, "general/records", "general_table", {}, records)
by = {r["_oid"].rsplit("-", 1)[-1]: r for r in records}
check("four values exported", len(records) == 4)
check("multi_reader published with label", by["0"]["reader_agreement"] == "multi_reader"
      and by["0"]["reader_agreement_label"] == "confirmed (two or more readers)")
check("single_reader published as 'unconfirmed (single reader)'", by["1"]["reader_agreement"] == "single_reader"
      and by["1"]["reader_agreement_label"] == "unconfirmed (single reader)")
check("single_reader_text_layer published", by["2"]["reader_agreement"] == "single_reader_text_layer"
      and by["2"]["reader_agreement_label"] == "single reader + PDF text layer")
check("no flag in the store -> field absent after clean()", "reader_agreement" not in E.clean(by["3"]))
check("second Claude reader keeps its id in readers", by["0"]["readers"] == ["claude-opus-5.5#2", "gemini-3.8-flash"])
check("first Claude reader unchanged", by["2"]["readers"] == ["claude-opus-5.5", "text_layer"])
check("CSV columns include reader_agreement and its label", "reader_agreement" in E.CSV_COLUMNS and "reader_agreement_label" in E.CSV_COLUMNS)
check("SCHEMA.md documents the field", "`reader_agreement`" in E.SCHEMA_MD and "unconfirmed (single reader)" in E.SCHEMA_MD)
print(f"{len(fails)} failures")
sys.exit(1 if fails else 0)
