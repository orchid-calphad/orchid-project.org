"""Table-context batches 1-15 (ORCHID-TC-20261027-v1): patch-file invariants and replay on the released records."""
import gzip, json, os, sys
from pathlib import Path
TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
import human_review_corrections as h
P = json.loads((TOOLS / "table_context_corrections_20261027_batches.json").read_text(encoding="utf-8"))
fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("" if ok else " -- %s" % (detail,)))
    if not ok: fails.append(name)
ids = [p["record_id"] for p in P]
check("one patch per record", len(ids) == len(set(ids)), len(ids) - len(set(ids)))
check("before-guards exclude gate-derived flags", not any("flags" in p["before"] for p in P))
check("patches change phases and structure only", all(set(p["after"]) <= {"phases", "structure"} for p in P),
      sorted({k for p in P for k in p["after"]} - {"phases", "structure"}))
check("every patch carries page evidence", all(p.get("evidence") and p.get("visual_check") for p in P))
lab = [ph for p in P for ph in (p["after"].get("phases") or [])]
check("structure fields name their printed source", all(ph.get("structure_basis") in (None, "printed_source_page") for ph in lab
                                                         if ph.get("pearson_symbol") or ph.get("prototype")))
lit = os.environ.get("ORCHID_LITERATURE", os.path.expanduser("~/projects/orchid-data/literature/orchid-literature.jsonl.gz"))
if os.path.exists(lit):
    recs = [json.loads(l) for l in gzip.open(lit, "rt", encoding="utf-8")]
    rep = h.apply(recs, P, evidence_id="ORCHID-TC-20261027-v1", gate="TC20261027", rule="source_adjudicated_table_context")
    check("released records already carry every patch (replay applies nothing)", rep["applied"] == [], len(rep["applied"]))
print("%d failed" % len(fails) if fails else "all passed")
sys.exit(1 if fails else 0)
