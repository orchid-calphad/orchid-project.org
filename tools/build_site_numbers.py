#!/usr/bin/env python3
"""Website numbers and citation blocks from the data (no hand-typed counts).

Every count on the hand-written website pages sits in an element marked `data-orchid="<key>"`, for example
`<b data-orchid="literature_values">77,540</b>`; this script recomputes each key from the data and rewrites the
element text. The release-wide "Cite" blocks (between `<!-- orchid-cite:begin ... -->` and `<!-- orchid-cite:end -->`)
of the pages in ALL_DATASET_PAGES list every dataset of tools/orchid_cite.py with its citation, the sources to cite as
well and the BibTeX; the page-specific blocks elsewhere are rebuilt the same way for the datasets they already cite
(found by their DOIs). Elements marked `data-orchid-name="<dataset>"` get the formal dataset name (ORCHID-Values, ...)
from the titles in tools/orchid_cite.py. The header, release status and footer of every page in site_chrome.PAGES
are rewritten from tools/site_chrome.py.

Reads (never writes) the orchid-data repository (literature, phase-diagrams, hulls, phases, assessed phase
assemblages, ESPEI bundle, bibliography) and this site's data/catalog.json. Deterministic; standard library only.

  python3 tools/build_site_numbers.py            # rewrite the marked numbers and cite blocks
  python3 tools/build_site_numbers.py --check    # exit 1 if a page is stale or names an unknown key
"""
from __future__ import annotations

import argparse
import glob
import html
import json
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402
import site_chrome  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGES = ("index.html", "dataset/index.html", "workflow/index.html", "about/index.html", "data-explorer/index.html",
         "phases/index.html", "phase-diagrams/index.html", "dft-parity/index.html")
ALL_DATASET_PAGES = ("index.html", "dataset/index.html", "workflow/index.html", "about/index.html")
WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
NAME = re.compile(r'(<(b|span|strong|h3)\b[^>]*\bdata-orchid-name="([a-z0-9-]+)"[^>]*>)(.*?)(</\2>)', re.S)
MARK = re.compile(r'(<(b|span|strong)\b[^>]*\bdata-orchid="([a-z0-9_]+)"[^>]*>)(.*?)(</\2>)', re.S)
CITE = re.compile(r"(<!-- orchid-cite:begin[^>]*-->\n)(.*?)(<!-- orchid-cite:end -->)", re.S)


def jload(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def n(x):
    return format(int(x), ",")


def numbers(data):
    """key -> display text, every value computed from the data."""
    out = {}
    cat = jload(os.path.join(SITE, "data", "catalog.json"))
    out["catalog_databases"] = n(len(cat))
    out["tdb_redistributable"] = n(sum(1 for c in cat if c.get("tier") == "bundled"))
    for layer, key in (("dft-mp", "hulls_mp"), ("dft-oqmd", "hulls_oqmd"), ("calphad-sgte", "hulls_sgte")):
        out[key] = n(len(glob.glob(os.path.join(data, "hulls", layer, "*", "*.json"))))
    t = jload(os.path.join(data, "literature", "index.json"))["totals"]
    out["literature_values"] = n(t["records"])
    out["literature_accepted"] = n(t["accepted"])
    out["literature_papers"] = n(t["papers_reporting"])
    out["literature_systems"] = n(t["systems"])
    out["sgte_tdb_files"] = n(sum(1 for c in cat if (c.get("local_path") or "").startswith("databases/sgte-binary/")))
    pdi = jload(os.path.join(data, "phase-diagrams", "index.json"))
    out["sgte_diagrams"] = n(sum(len(v) for v in pdi["systems"].values()))
    reg = jload(os.path.join(data, "phases", "index.json"))
    rc = reg.get("counts") or {}
    out["registry_phases"] = n(rc.get("phases", 0))
    out["registry_systems"] = n(rc.get("systems", 0))
    asm = jload(os.path.join(data, "assessed-phase-assemblages", "summary.json"))
    out["assemblage_points"] = n(asm["points"])
    out["assemblage_diagrams"] = n(asm["diagrams_with_points"])
    es = jload(os.path.join(data, "espei", "index.json"))["totals"]
    pc = (es.get("phase_constraints") or {}).get("folders") or {}
    pc_sets, pc_pts = sum(v["datasets"] for v in pc.values()), sum(v["points"] for v in pc.values())
    # core: literature- and DFT-derived datasets; constraint: optional phase-constraint datasets (assessed handbook
    # assemblages, derived invariant brackets and liquid-stability points); all: both, as released in espei/index.json
    out["espei_datasets"] = n(es["datasets"])
    out["espei_values"] = n(es["values"])
    out["espei_constraint_datasets"] = n(pc_sets)
    out["espei_constraint_points"] = n(pc_pts)
    out["espei_datasets_all"] = n(es["datasets"] + pc_sets)
    out["espei_values_all"] = n(es["values"] + pc_pts)
    bc = jload(os.path.join(data, "bibliography", "index.json"))["counts"]
    out["bibliography_papers"] = n(bc["papers"])
    out["bibliography_edges"] = n(bc["edges"])
    out["datasets_word"] = WORDS[len(orchid_cite.NAMES) - 1]
    return out


def formal_name(name):
    """Formal dataset name: the part of the orchid_cite title before the colon (ORCHID-Values, ORCHID ESPEI datasets)."""
    return orchid_cite.DATASETS[name][0].split(":", 1)[0]


def cited_names(old, data):
    """Datasets a page-specific block already cites, found by DOI (placeholder or written back)."""
    by_doi = {v: k for k, v in orchid_cite.dois(data).items()}
    return [nm for nm in orchid_cite.NAMES if nm in {by_doi.get(d) for d in re.findall(r"doi\.org/(10\.[^\"<\s]+)", old)}]


def cite_block(old, data, names=None):
    """A cite block: the given datasets (default every dataset), then the umbrella; heading and intro kept."""
    names = [x for x in (names or orchid_cite.NAMES) if x != orchid_cite.UMBRELLA] + [orchid_cite.UMBRELLA]
    head = old.split("<ul>", 1)[0]
    head = re.sub(r"\b(one|two|three|four|five|six|seven|eight|nine|ten) (citable )?datasets",
                  lambda m: WORDS[len(names) - 1] + " " + (m.group(2) or "") + "datasets", head)
    lis = []
    for nm in names:
        e = orchid_cite.entry(nm, data)
        text = e["text"].rsplit(" https://doi.org/", 1)[0]
        lis.append('<li>%s <a href="%s">doi:%s</a></li>' % (html.escape(text, quote=False), e["url"], e["doi"]))
    also = []
    for nm in names:
        for k in orchid_cite.DATASETS[nm][2]:
            if orchid_cite.THIRD_PARTY[k] not in also:
                also.append(orchid_cite.THIRD_PARTY[k])
    return (head + "<ul>\n" + "\n".join(lis) + "\n</ul>\n" +
            '<p class="muted">Also cite %s.</p>\n' % html.escape("; ".join(also), quote=False) +
            "<pre>%s</pre>\n</section>\n" % html.escape(orchid_cite.bibtex_block(names[:-1], data)))


def render(path, text, nums, data):
    unknown = []

    def sub(m):
        key = m.group(3)
        if key not in nums:
            unknown.append(key)
            return m.group(0)
        return m.group(1) + nums[key] + m.group(5)
    def name_sub(m):
        if m.group(3) not in orchid_cite.DATASETS:
            unknown.append(m.group(3))
            return m.group(0)
        return m.group(1) + html.escape(formal_name(m.group(3)), quote=False) + m.group(5)
    new = MARK.sub(sub, text)
    new = NAME.sub(name_sub, new)
    if path in ALL_DATASET_PAGES:
        new = CITE.sub(lambda m: m.group(1) + cite_block(m.group(2), data) + m.group(3), new)
    else:
        new = CITE.sub(lambda m: m.group(1) + cite_block(m.group(2), data, cited_names(m.group(2), data)) + m.group(3), new)
    return new, unknown


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data repository (read only)")
    ap.add_argument("--check", action="store_true", help="exit 1 if a page is stale")
    a = ap.parse_args(argv)
    nums = numbers(a.data)
    bad, wrote, used = [], [], set()
    for rel in PAGES + tuple(p for p in site_chrome.PAGES if p not in PAGES):
        p = os.path.join(SITE, rel)
        text = open(p, encoding="utf-8").read()
        if rel in PAGES:
            used.update(m.group(3) for m in MARK.finditer(text))
            new, unknown = render(rel, text, nums, a.data)
            bad += ["%s: unknown key %s" % (rel, k) for k in unknown]
        else:
            new = text
        if rel in site_chrome.PAGES:
            new, missing = site_chrome.apply(rel, new)
            bad += ["%s: chrome marker missing: %s" % (rel, k) for k in missing]
        if new != text:
            if a.check:
                bad.append("stale: " + rel)
            else:
                with open(p, "w", encoding="utf-8") as fh:
                    fh.write(new)
                wrote.append(rel)
    if a.check:
        for b in bad:
            print(b)
        print("FAIL: %d problem(s)" % len(bad) if bad else "OK: %d pages current (%d keys used; chrome on %d pages)" % (
            len(PAGES), len(used), len(site_chrome.PAGES)))
        return 1 if bad else 0
    for b in bad:
        print(b)
    print("wrote %s" % (", ".join(wrote) or "nothing (current)"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
