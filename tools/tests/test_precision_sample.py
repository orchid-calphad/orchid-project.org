#!/usr/bin/env python3
"""Frozen precision sample (tools/build_precision_sample.py, tools/score_precision_sample.py).

Fixture: a small literature export is drawn into a temp directory; then the export moves on. Checks: `draw` refuses a
directory that holds a sample (and the kit's own directory); `--check` passes while the frozen files are untouched,
fails when sample.json / index.html / a packet changes or a sampled id disappears, and only REPORTS (no failure)
values changed since sampling; ids merged by the de-duplication count as released; `freeze` refuses to rewrite a
freeze record; the scorer refuses an unfrozen or altered sample and reports per item whether the release still carries
the reviewed value; zenodo_release.py only ever runs the builder with --check. Also, when present, the distributed
kit (review/precision-sample) is frozen at 5c3d02e5... and verifies against its population. Temp files only.
usage: python3 tools/tests/test_precision_sample.py"""
import contextlib, gzip, hashlib, io, json, os, re, shutil, sys, tempfile
from pathlib import Path
sys.dont_write_bytecode = True
TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))
import build_precision_sample as B  # noqa: E402
import score_precision_sample as S  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            return fn(*a, **k)
        except SystemExit as e:
            return ("exit", str(e))


TMP = Path(tempfile.mkdtemp(prefix="orchid-ps-"))
lit = TMP / "lit"
lit.mkdir()


def rec(i, ev="experimental", kind="general_table", value=None):
    return {"id": "ORC-%05d" % i, "system": "AL-NI", "record_kind": kind, "evidence_class": ev, "property": "HM_FORM",
            "value": -40.0 - i if value is None else value, "unit": "kJ/mol-atom", "paper": "doi:10.1000/p%d" % (i % 7),
            "phases": [{"label": "AlNi"}], "flags": {"extraction_accepted": True}}


def write_release(recs, papers):
    (lit / "orchid-literature.jsonl.gz").write_bytes(gzip.compress("".join(json.dumps(r) + "\n" for r in recs).encode(), mtime=0))
    (lit / "papers.json").write_text(json.dumps(papers), encoding="utf-8")
    (lit / "index.json").write_text(json.dumps({"release": "fixture"}), encoding="utf-8")


recs = [rec(i, ev=("experimental", "model", "dft")[i % 3], kind=B.KINDS[i % 3]) for i in range(400)]
papers = {"doi:10.1000/p%d" % j: {"doi": "10.1000/p%d" % j, "citation": "P%d" % j} for j in range(7)}
write_release(recs, papers)
kit = TMP / "kit"


def run(*args):
    return quiet(B.main, list(args) + ["--input", str(lit / "orchid-literature.jsonl.gz"), "--papers", str(lit / "papers.json"),
                                       "--out", str(kit)])


check("draw into a new directory", run("draw", "-n", "30") == 0)
check("draw writes the freeze record", (kit / "frozen.json").exists())
check("draw refuses an existing sample", run("draw", "-n", "30") == 1)
check("draw refuses the distributed kit directory",
      quiet(B.main, ["draw", "--out", str(B.OUTDIR)]) == 1)
(kit / "packets").mkdir()
(kit / "packets" / "p1.html").write_text("packet", encoding="utf-8")
check("check flags a packet missing from the freeze record", run("--check") == 1)
fz = json.loads((kit / "frozen.json").read_text(encoding="utf-8"))
check("freeze refuses to rewrite a different freeze record", run("freeze") == 1)
(kit / "packets" / "p1.html").unlink()
check("check passes on the untouched kit", run("--check") == 0)
check("default command is check (verify only)", run() == 0)
sampled = [it["id"] for it in fz["items"]]
check("freeze record lists every sampled id with reviewed values and hashes",
      len(sampled) == len(json.loads((kit / "sample.json").read_text(encoding="utf-8"))["items"]) and all(it["reviewed_sha256"] == B.view_sha(it["reviewed"]) for it in fz["items"]))
check("population snapshot recorded", fz["population"]["sha256"] == hashlib.sha256(
    (lit / "orchid-literature.jsonl.gz").read_bytes()).hexdigest() and fz["population"]["population_accepted"] == 400)

# the release moves on: +values, one sampled value changed, one sampled id merged into another, citation DOI changed
before = (kit / "sample.json").read_bytes()
by_id = {r["id"]: r for r in recs}
chg, mrg, cit = sampled[0], sampled[1], sampled[2]
by_id[chg]["value"] = 999.0
survivor = dict(by_id.pop(mrg), id="ORC-SURVIVOR", merged_ids=[mrg])
by_id["ORC-SURVIVOR"] = survivor
pk = by_id[cit]["paper"]
papers2 = dict(papers)
papers2[pk] = dict(papers[pk], doi=None, doi_status="dropped_mismatch")
new = list(by_id.values()) + [rec(1000 + i) for i in range(50)]
write_release(new, papers2)
drift = TMP / "drift.json"
check("check still passes after the export changed (sample is frozen)", run("--check", "--drift-json", str(drift)) == 0)
check("sample.json untouched by check", (kit / "sample.json").read_bytes() == before)
d = {x["id"]: x for x in json.loads(drift.read_text(encoding="utf-8"))["items"]}
check("changed value reported (value question)", d[chg]["status"] == "changed" and d[chg]["questions"] == ["value"])
check("merged id counts as released", d[mrg].get("merged_into") == "ORC-SURVIVOR" and d[mrg]["status"] != "missing")
others_on_paper = [i for i in sampled if by_id.get(i, {}).get("paper") == pk]
check("dropped paper DOI reported (citation question)", all(d[i]["status"] == "changed" and "citation" in d[i]["questions"]
                                                             for i in others_on_paper) and others_on_paper)
# a sampled id disappears -> failure
gone = sampled[3]
write_release([r for r in new if r["id"] != gone], papers2)
check("check fails when a sampled id is no longer released", run("--check") == 1)
write_release(new, papers2)
# tampering with frozen files -> failure
(kit / "sample.json").write_bytes(before + b" ")
check("check fails when sample.json changed", run("--check") == 1)
(kit / "sample.json").write_bytes(before)
html = (kit / "index.html").read_bytes()
(kit / "index.html").write_bytes(html.replace(b"Precision sample", b"Precision  sample", 1))
check("check fails when index.html changed", run("--check") == 1)
(kit / "index.html").write_bytes(html)
check("check passes again once restored", run("--check") == 0)

# scorer: frozen sample only, per-item release status
s = json.loads(before)
sha = hashlib.sha256(before).hexdigest()
ans = [{"id": it["record"]["id"], "n": it["n"], "stratum": it["stratum"], "value": "yes", "unit": "yes",
        "evidence_class": "yes", "phases": "yes", "citation": "yes", "note": "", "updated": None} for it in s["items"][:5]]
res = TMP / "res.json"
res.write_text(json.dumps({"schema": "orchid-precision-review/1", "sample_sha256": sha, "reviewer": "fixture", "answers": ans}))
out_json = TMP / "scores.json"
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rc = S.main([str(res), "--sample", str(kit / "sample.json"), "--input", str(lit / "orchid-literature.jsonl.gz"),
                 "--papers", str(lit / "papers.json"), "--json", str(out_json)])
txt = buf.getvalue()
sc = json.loads(out_json.read_text())
check("scorer runs on the frozen sample", rc == 0 and "Release check" in txt)
check("scorer reports per item whether the release carries the reviewed value",
      sc["release_check"][chg]["status"] == "changed" and sc["release_check"][mrg]["merged_into"] == "ORC-SURVIVOR"
      and set(sc["release_check"]) == {a["id"] for a in ans})
fz_path = kit / "frozen.json"
fz_bytes = fz_path.read_bytes()
fz_path.unlink()
check("scorer refuses an unfrozen sample", quiet(S.main, [str(res), "--sample", str(kit / "sample.json"), "--input",
                                                          str(lit / "orchid-literature.jsonl.gz"), "--papers", str(lit / "papers.json")])[0] == "exit")
fz_path.write_bytes(fz_bytes)

# the release tool never regenerates the sample
src = (TOOLS / "zenodo_release.py").read_text(encoding="utf-8")
calls = re.findall(r"\[[^\[\]]*build_precision_sample\.py[^\[\]]*\]", src)
check("zenodo_release.py calls build_precision_sample.py only with --check", calls and all('"--check"' in c for c in calls))
check("zenodo_release.py text replacement skips review/", '"review/"' in src)

# the distributed kit
real = B.OUTDIR / "frozen.json"
if real.exists():
    fz = json.loads(real.read_text(encoding="utf-8"))
    check("distributed kit frozen at 5c3d02e5...", fz["sample_sha256"].startswith("5c3d02e5")
          and hashlib.sha256((B.OUTDIR / "sample.json").read_bytes()).hexdigest() == fz["sample_sha256"])
    check("distributed kit: re-draw from its population reproduced it", fz["population"]["draw_reproduced"] is True)
    bad, _ = B.verify(B.OUTDIR, B.INPUT, B.PAPERS)
    check("distributed kit verifies against the current export", bad == [])
shutil.rmtree(TMP, ignore_errors=True)
print("%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
