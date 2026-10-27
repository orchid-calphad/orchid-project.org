#!/usr/bin/env python3
"""zenodo_release.py writeback safety (temporary directories only; never touches the real repositories).

1. sync_literature never deletes: a file present only in <data>/literature is moved to
   <superseded>/<stamp>-zenodo-writeback/literature/<path>, listed in that folder's README.md; changed files are copied.
2. writeback refuses a partial DOI map (all ten records with ten distinct real DOIs are required), and after the
   substitution it searches the whole candidate tree for any remaining `zenodo.PENDING-` placeholder, including names
   that are not in the map and placeholders inside .gz files, and fails if one remains.
usage: python3 tools/tests/test_zenodo_writeback.py"""
import argparse, gzip, json, os, shutil, sys, tempfile
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import orchid_cite as oc  # noqa: E402
import zenodo_release as Z  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def write(p, text, gz=False):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with (gzip.open(p, "wt", encoding="utf-8") if gz else open(p, "w", encoding="utf-8")) as fh:
        fh.write(text)


tmp = tempfile.mkdtemp(prefix="orchid-writeback-test-")
try:
    # ---------------------------------------------------------------------------------------- 1. sync --
    src, data, sup = (os.path.join(tmp, x) for x in ("site-literature", "data", "superseded"))
    write(os.path.join(src, "index.json"), "new index")
    write(os.path.join(src, "systems", "AL-NI.json"), "same")
    write(os.path.join(data, "literature", "index.json"), "old index")
    write(os.path.join(data, "literature", "systems", "AL-NI.json"), "same")
    write(os.path.join(data, "literature", "systems", "OLD-ONLY.json"), "only in destination")
    write(os.path.join(data, "literature", "stale.md"), "stale note")
    os.makedirs(os.path.join(sup, "2026-01-01-zenodo-writeback", "literature"))
    write(os.path.join(sup, "2026-01-01-zenodo-writeback", "literature", "stale.md"), "earlier superseded copy")
    moved = Z.sync_literature(data, src=src, superseded_root=sup, stamp="2026-01-01")
    box = os.path.join(sup, "2026-01-01-zenodo-writeback")
    check("changed file copied", open(os.path.join(data, "literature", "index.json")).read() == "new index")
    check("destination-only file no longer in the data tree", not os.path.exists(os.path.join(data, "literature", "systems", "OLD-ONLY.json")))
    check("destination-only file moved, content intact",
          open(os.path.join(box, "literature", "systems", "OLD-ONLY.json")).read() == "only in destination")
    check("earlier superseded copy not overwritten",
          open(os.path.join(box, "literature", "stale.md")).read() == "earlier superseded copy"
          and open(os.path.join(box, "literature", "stale.md.1")).read() == "stale note")
    readme = open(os.path.join(box, "README.md")).read()
    check("README lists every moved file", "literature/systems/OLD-ONLY.json" in readme and "literature/stale.md.1" in readme)
    check("two files reported moved", len(moved) == 2)
    check("nothing left to move on a second sync", Z.sync_literature(data, src=src, superseded_root=sup, stamp="2026-01-01") == [])

    # ------------------------------------------------------------------------------------ 2. writeback --
    names = list(Z.ALL_RECORDS)
    check("ten records", len(names) == 10)
    real = {n: "10.5281/zenodo.%d" % (9000000 + i) for i, n in enumerate(names)}

    def tree():
        root = os.path.join(tmp, "wb")
        shutil.rmtree(root, ignore_errors=True)
        d, s = os.path.join(root, "data"), os.path.join(root, "site")
        write(os.path.join(d, "zenodo", "dois.json"), json.dumps({"dois": {n: oc.placeholder(n) for n in names}}))
        write(os.path.join(d, "README.md"), " ".join("https://doi.org/" + oc.placeholder(n) for n in names))
        write(os.path.join(s, "index.html"), "<a>%s</a>" % oc.placeholder("literature"))
        return d, s

    def ns(d, s, dois, **kw):
        return argparse.Namespace(data=d, site=s, examples=os.path.join(tmp, "no-examples"), examples_ref=None,
                                  doi=dois, dry_run=kw.get("dry_run", False), no_rebuild=True, parquet_python=sys.executable)

    d, s = tree()
    partial = ["%s=%s" % (n, real[n]) for n in names[:9]]
    try:
        Z.cmd_writeback(ns(d, s, partial))
        refused = False
    except SystemExit as e:
        refused = "missing" in str(e) and names[9] in str(e)
    check("partial DOI map refused (names the missing record)", refused)
    check("partial map wrote nothing", oc.placeholder(names[0]) in open(os.path.join(d, "README.md")).read())
    try:
        Z.cmd_writeback(ns(d, s, partial, dry_run=True))
        refused = False
    except SystemExit:
        refused = True
    check("partial DOI map refused in a dry run too", refused)
    dup = ["%s=%s" % (n, real[names[0]]) for n in names]
    try:
        Z.cmd_writeback(ns(d, s, dup))
        refused = False
    except SystemExit:
        refused = True
    check("non-distinct DOIs refused", refused)
    try:
        Z.cmd_writeback(ns(d, s, ["%s=%s" % (n, real[n]) for n in names[:9]] + ["%s=%s" % (names[9], oc.placeholder(names[9]))]))
        refused = False
    except SystemExit:
        refused = True
    check("placeholder given as a DOI refused", refused)

    full = ["%s=%s" % (n, real[n]) for n in names]
    check("complete map, clean tree: writeback passes", Z.cmd_writeback(ns(d, s, full)) == 0)
    check("all placeholders replaced", not Z.find_placeholders(d) and not Z.find_placeholders(s))
    check("dois.json holds the ten DOIs", json.load(open(os.path.join(d, "zenodo", "dois.json")))["dois"] == real)

    d, s = tree()
    write(os.path.join(d, "notes.md"), "stray 10.5281/zenodo.PENDING-unknown-record placeholder")
    write(os.path.join(d, "literature", "dump.jsonl.gz"), '{"doi": "10.5281/zenodo.PENDING-literature"}\n', gz=True)
    write(os.path.join(s, "tools", "orchid_cite_copy.py"), oc.PLACEHOLDER_PREFIX)   # tools/ is code: skipped
    rc = Z.cmd_writeback(ns(d, s, full))
    left = Z.find_placeholders(d)
    check("leftover placeholder fails writeback", rc == 1)
    check("leftover outside the map found", "notes.md" in left)
    check("leftover inside a .gz found", "literature/dump.jsonl.gz" in left)
    check("tools/ code not reported", not Z.find_placeholders(s))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
