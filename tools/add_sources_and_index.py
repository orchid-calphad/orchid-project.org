#!/usr/bin/env python3
"""Attach original-assessment citation (TDB header) + DOI (SGTE-corpus tdb_pdf_links.json) to each diagram JSON and
rebuild data/phase-diagrams/index.json (2026-10-22).

A header citation whose matched DOI is a mismatch (tools/doi_corrections.json, context tdb_header) keeps the printed
citation with `doi: null`, `doi_status: "dropped_mismatch"` and `doi_correction: <entry id>` (2026-10-23).

  python3 tools/add_sources_and_index.py [ROOT ...]     # default: this site's data/ copy and the orchid-data repository
  python3 tools/add_sources_and_index.py --check        # exit 1 if a diagram file or index.json would change
Run tools/add_cite_as.py afterwards (it adds `cite_as` to index.json)."""
import collections, glob, json, os, re, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import doi_corrections  # noqa: E402
H = os.path.expanduser("~"); SITE = f"{H}/projects/orchid-project.org"; T = f"{H}/projects/open-tdb/databases/sgte-binary"
DATA = os.environ.get("ORCHID_DATA") or os.path.join(os.path.dirname(SITE), "orchid-data")
CORR = doi_corrections.load()
L = {x["tdb"]: x for x in json.load(open(f"{H}/projects/open-tdb/build/sgte-corpus/tdb_pdf_links.json"))["links"]}
def header(tdb):
    lines = []
    for l in open(f"{T}/{tdb}", errors="replace"):
        if not l.startswith("$"): break
        lines.append(l[1:].strip())
    m = re.search(r"Database for .*? from (.*?)\s+Dataset created", " ".join(lines))
    c = (m.group(1) if m else "").strip()
    return c.split("  ")[-1].strip() if "  " in c else c
def build(root):
    """{path: text} for every diagram file with status ok and the index (without cite_as)."""
    out = {}; sysd = collections.defaultdict(list); n = ok = nd = nc = 0
    for f in sorted(glob.glob(f"{root}/phase-diagrams/sgte/*.json")):
        d = json.load(open(f)); n += 1
        if d.get("status") != "ok": continue
        x = L.get(d["tdb"], {}); src = {"citation": header(d["tdb"]) or x.get("header_citation", ""), "doi": x.get("doi") or None}
        e = CORR.tdb(d["tdb"], src["doi"]) if src["doi"] else None
        if e: CORR.mark(src, e); nc += 1
        if src["doi"]: src["url"] = "https://doi.org/" + src["doi"]; nd += 1
        d["original_source"] = src; out[f] = json.dumps(d, separators=(",", ":")); ok += 1
        row = {"file": os.path.basename(f), "tdb": d["tdb"], "sha256": d["sha256"], "grid": "grid" in d,
               "citation": src["citation"], "doi": src["doi"]}
        if e: row.update(doi_status=src["doi_status"], doi_correction=src["doi_correction"])
        sysd["-".join(d["components"])].append(row)
    out[f"{root}/phase-diagrams/index.json"] = json.dumps(
        {"source": "SGTE binary collection (Hallstedt 2025), CC BY 4.0", "engine": "pycalphad 0.11.2 BinaryStrategy + equilibrium grid",
         "systems": dict(sorted(sysd.items()))}, indent=1, ensure_ascii=False)
    return out, f"{root}: {n} files, {ok} ok, {nd} with DOI, {nc} DOI dropped as mismatch, {len(sysd)} systems"
def main(argv):
    check = "--check" in argv
    roots = [a for a in argv if not a.startswith("--")] or [os.path.join(SITE, "data"), DATA]
    bad = []
    for root in roots:
        out, msg = build(root)
        for p, text in out.items():
            cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
            if p.endswith("index.json") and cur is not None:   # compare without the cite_as added by add_cite_as.py
                cur = json.dumps({k: v for k, v in json.loads(cur).items() if k != "cite_as"}, indent=1, ensure_ascii=False)
            if cur != text:
                if check: bad.append(p); continue
                if p.endswith("index.json") and os.path.exists(p):   # keep an existing cite_as (add_cite_as.py refreshes it)
                    old = json.load(open(p, encoding="utf-8"))
                    if "cite_as" in old:
                        dd = json.loads(text); dd["cite_as"] = old["cite_as"]; text = json.dumps(dd, indent=1, ensure_ascii=False)
                open(p, "w", encoding="utf-8").write(text)
        print(("check " if check else "index: ") + msg)
    if bad:
        print("stale: " + ", ".join(bad[:10]) + (" ... (%d)" % len(bad) if len(bad) > 10 else ""), file=sys.stderr)
        return 1
    return 0
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
