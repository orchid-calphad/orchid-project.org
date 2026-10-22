#!/usr/bin/env python3
"""Attach original-assessment citation (TDB header) + DOI (SGTE-corpus tdb_pdf_links.json) to each diagram JSON and
rebuild data/phase-diagrams/index.json (2026-10-22)."""
import collections, glob, json, os, re
H = os.path.expanduser("~"); SITE = f"{H}/projects/orchid-project.org"; T = f"{H}/projects/open-tdb/databases/sgte-binary"
L = {x["tdb"]: x for x in json.load(open(f"{H}/projects/open-tdb/build/sgte-corpus/tdb_pdf_links.json"))["links"]}
def header(tdb):
    lines = []
    for l in open(f"{T}/{tdb}", errors="replace"):
        if not l.startswith("$"): break
        lines.append(l[1:].strip())
    m = re.search(r"Database for .*? from (.*?)\s+Dataset created", " ".join(lines))
    c = (m.group(1) if m else "").strip()
    return c.split("  ")[-1].strip() if "  " in c else c
sysd = collections.defaultdict(list); n = ok = nd = 0
for f in sorted(glob.glob(f"{SITE}/data/phase-diagrams/sgte/*.json")):
    d = json.load(open(f)); n += 1
    if d.get("status") != "ok": continue
    x = L.get(d["tdb"], {}); src = {"citation": header(d["tdb"]) or x.get("header_citation", ""), "doi": x.get("doi") or None}
    if src["doi"]: src["url"] = "https://doi.org/" + src["doi"]; nd += 1
    d["original_source"] = src; json.dump(d, open(f, "w"), separators=(",", ":")); ok += 1
    sysd["-".join(d["components"])].append({"file": os.path.basename(f), "tdb": d["tdb"], "sha256": d["sha256"],
        "grid": "grid" in d, "citation": src["citation"], "doi": src["doi"]})
json.dump({"source": "SGTE binary collection (Hallstedt 2025), CC BY 4.0", "engine": "pycalphad 0.11.2 BinaryStrategy + equilibrium grid",
           "systems": dict(sorted(sysd.items()))}, open(f"{SITE}/data/phase-diagrams/index.json", "w"), indent=1, ensure_ascii=False)
print(f"index: {n} files, {ok} ok, {nd} with DOI, {len(sysd)} systems")
