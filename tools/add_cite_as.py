#!/usr/bin/env python3
"""Add the `cite_as` object and per-folder CITATION.cff files to the data that no other builder here writes.

The calculated phase-diagram index (written by tools/add_sources_and_index.py) and the hull catalogs (copied from
open-tdb, where tools/build_binary_hulls.py writes them) get a top-level `cite_as` (tools/orchid_cite.py); the
rest of each file is re-serialised byte-for-byte in its original style (the script refuses a file whose
style it cannot reproduce). It also writes

  phase-diagrams/CITATION.cff   hulls/CITATION.cff

Literature files (index.json, systems/*.json, CITATION.txt, CITATION.cff) are written by
tools/export_literature_data.py, dft_parity.json by tools/build_dft_parity.py, and Parquet by tools/build_parquet.py.

Applied to every root given (default: this site's data/ copy and the orchid-data repository) so the two stay
identical. DOIs come from <orchid-data>/zenodo/dois.json. Idempotent; standard library only; no network.

  python3 tools/add_cite_as.py            # write
  python3 tools/add_cite_as.py --check    # exit 1 if any file is missing or stale
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def style_pd(d):  # tools/add_sources_and_index.py: indent 1, insertion order, no trailing newline
    return json.dumps(d, indent=1, ensure_ascii=False)


def style_hull(d):  # open-tdb build_binary_hulls.py: indent 2, sorted keys, trailing newline
    return json.dumps(d, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


TARGETS = [  # (relative path, dataset, serialiser)
    ("phase-diagrams/index.json", "phase-diagrams", style_pd),
    ("hulls/dft-mp/catalog.json", "hulls", style_hull),
    ("hulls/dft-oqmd/catalog.json", "hulls", style_hull),
    ("hulls/calphad-sgte/catalog.json", "hulls", style_hull),
]
CFF = [("phase-diagrams/CITATION.cff", "phase-diagrams", " (the SGTE binary collection, Hallstedt 2025, and the "
        "original assessment named in each TDB file)"),
       ("hulls/CITATION.cff", "hulls", " (Materials Project, OQMD and the SGTE binary collection)")]


def expected(root, dois_root):
    out = {}
    for rel, ds, style in TARGETS:
        p = os.path.join(root, rel)
        with open(p, encoding="utf-8") as fh:
            raw = fh.read()
        d = json.loads(raw)
        base = {k: v for k, v in d.items() if k != "cite_as"}
        if "cite_as" not in d and style(base) != raw:
            raise SystemExit("refusing %s: cannot reproduce its JSON style byte-for-byte" % p)
        if style is style_pd:  # keep the original key order; cite_as first so it is visible at the top
            new = {"cite_as": orchid_cite.cite_as(ds, dois_root), **base}
        else:
            new = dict(base, cite_as=orchid_cite.cite_as(ds, dois_root))
        out[rel] = style(new)
    for rel, ds, extra in CFF:
        out[rel] = orchid_cite.cff(ds, dois_root, extra)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", action="append", help="data root holding phase-diagrams/ and hulls/ (repeatable; "
                    "default: <site>/data and the orchid-data repository)")
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data repository (zenodo/dois.json)")
    ap.add_argument("--check", action="store_true", help="exit 1 if any file is missing or stale")
    a = ap.parse_args(argv)
    roots = a.root or [os.path.join(SITE, "data"), a.data]
    bad, wrote = [], 0
    for root in roots:
        if not os.path.isdir(os.path.join(root, "phase-diagrams")):
            print("skip (no phase-diagrams/): " + root, file=sys.stderr)
            continue
        for rel, text in expected(root, a.data).items():
            p = os.path.join(root, rel)
            cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
            if cur == text:
                continue
            if a.check:
                bad.append(p)
            else:
                with open(p + ".tmp", "w", encoding="utf-8") as fh:
                    fh.write(text)
                os.replace(p + ".tmp", p)
                wrote += 1
    if a.check:
        if bad:
            print("stale:\n  " + "\n  ".join(bad), file=sys.stderr)
            return 1
        print("cite_as current in %d root(s)" % len(roots))
        return 0
    print("wrote %d file(s) in %d root(s)" % (wrote, len(roots)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
