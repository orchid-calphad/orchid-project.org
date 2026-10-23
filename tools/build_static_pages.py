#!/usr/bin/env python3
"""Build static, JavaScript-free pages of the ORCHID data (readable by people, search engines and AI crawlers).

Reads (never modifies) the data repository (default: ../orchid-data next to this site repository):

  literature/index.json, literature/systems/<SYS>[.partN].json, literature/papers.json, literature/dft_parity.json
  phase-diagrams/index.json, phase-diagrams/sgte/*.json, hulls/{dft-mp,dft-oqmd,calphad-sgte}/

and writes, into the data repository (served at https://data.orchid-project.org/):

  systems/<SYS>.html|.md|.csv   one page per chemical system with literature data and/or a calculated diagram
  systems/index.html|.md        list of systems
  compounds/<FORMULA>.html|.md  formation enthalpies per compound (literature, Materials Project, OQMD)
  compounds/index.html|.md      list of compounds
  papers/index.html|.md         reporting papers (citation, DOI, systems, number of values)
  index.html|.md                what the data site is, folder guide, licence, citation
  sitemap.xml, robots.txt, llms.txt, llms-full.txt

and, into this site repository: sitemap.xml (sitemap index), sitemap-pages.xml, robots.txt.

Published content: derived data, citations / DOIs, and the one verbatim sentence per value already in the
literature release (fair use). No page images, crops or full text. All text is HTML-escaped. Pages carry
schema.org Dataset JSON-LD and a Markdown twin. Deterministic: same inputs and --date give identical bytes.
Standard library only; no network.

  python3 tools/build_static_pages.py                    # write
  python3 tools/build_static_pages.py --check            # exit 1 if any output is missing, stale or orphaned
  python3 tools/build_static_pages.py --date 2026-10-22  # lastmod date written to the sitemaps
"""
from __future__ import annotations

import argparse
import csv
import html
import importlib.util
import io
import json
import os
import posixpath
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction

sys.dont_write_bytecode = True  # importing build_dft_parity must not leave __pycache__ behind
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA = os.path.join(os.path.dirname(SITE), "orchid-data")
DATA_URL = "https://data.orchid-project.org/"
SITE_URL = "https://orchid-project.org/"
LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"
CITE_DATA = DEFAULT_DATA  # repository holding zenodo/dois.json (set from --data in main)
CAP_VALUES = 1500
QUOTE_MAX = 400
DATASET_NAME = "ORCHID data: literature thermodynamic data, calculated phase diagrams and DFT hulls"
OWN_DIRS = ("systems", "compounds", "papers")
ROOT_FILES = ("index.html", "index.md", "sitemap.xml", "robots.txt", "llms.txt", "llms-full.txt")
AI_BOTS = ("GPTBot", "ChatGPT-User", "OAI-SearchBot", "ClaudeBot", "Claude-Web", "Claude-User", "Claude-SearchBot",
           "anthropic-ai", "PerplexityBot", "Perplexity-User", "Google-Extended", "CCBot", "Applebot-Extended",
           "Bytespider", "Amazonbot", "meta-externalagent", "cohere-ai", "DuckAssistBot")
SITE_PAGES = ("", "phase-diagrams/", "data-explorer/", "dft-parity/", "dataset/", "about/", "workflow/", "catalog/",
              "catalog/binary-thermo.html", "catalog/mixing-energy.html", "catalog/model-audit/",
              "catalog/compatibility/unary-compatibility-network.html", "best-practices.html")

EVL = {"experimental": "experimental", "model": "model", "dft": "DFT", "unspecified": "unspecified"}
EV_ORDER = ("experimental", "model", "dft", "unspecified")
PROPN = {"INVARIANT_TEMPERATURE": "Invariant reaction", "TRANSITION_TEMPERATURE": "Phase transition",
         "PHASE_COMPOSITION": "Phase composition", "PHASE_ASSEMBLAGE": "Phase assemblage",
         "HM_FORM": "Enthalpy of formation", "SM_FORM": "Entropy of formation", "GM_FORM": "Gibbs energy of formation",
         "HM_MIX": "Enthalpy of mixing", "GM_MIX": "Gibbs energy of mixing", "SM_MIX": "Entropy of mixing",
         "CPM": "Heat capacity", "SM": "Entropy", "ENTROPY": "Entropy", "HEAT_CONTENT": "Heat content",
         "GIBBS_ENERGY": "Gibbs energy", "PARTIAL_ENTHALPY": "Partial enthalpy", "EMF": "EMF", "ACTIVITY": "Activity",
         "ACTIVITY_COEFFICIENT": "Activity coefficient", "CHEMICAL_POTENTIAL": "Chemical potential",
         "VAPOR_PRESSURE": "Vapour pressure", "MOLAR_VOLUME": "Molar volume",
         "TRANSITION_ENTHALPY": "Transition enthalpy", "LATTICE_A": "Lattice a", "LATTICE_B": "Lattice b",
         "LATTICE_C": "Lattice c", "OTHER": "Other"}
GROUP_ORDER = ("invariant reaction", "phase transition", "phase equilibrium", "formation", "mixing / partial",
               "heat capacity / entropy / enthalpy", "activity / EMF / Gibbs energy", "crystal structure",
               "model expression", "other property")

CSS = """:root{--bg:#fbfaf7;--fg:#1d1b20;--muted:#5e5a66;--line:#e2dfe6;--card:#fff;--accent:#7b3fa0;--row:#f6f4f8;
--exp:#1d7448;--exp-bg:#e3f3ea;--mod:#2c58a8;--mod-bg:#e5ecf9;--dft:#a3520a;--dft-bg:#fbeedf;--uns:#5f5a66;--uns-bg:#ecebee}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#141218;--fg:#ece8f1;--muted:#a69fb0;--line:#2e2a35;
--card:#1c1922;--accent:#c79ae6;--row:#221e28;--exp:#74d3a2;--exp-bg:#183326;--mod:#9dbcf3;--mod-bg:#1b2840;--dft:#f2ad6b;
--dft-bg:#3a2612;--uns:#bdb7c6;--uns-bg:#2a2730}}
:root[data-theme="dark"]{--bg:#141218;--fg:#ece8f1;--muted:#a69fb0;--line:#2e2a35;--card:#1c1922;--accent:#c79ae6;--row:#221e28;
--exp:#74d3a2;--exp-bg:#183326;--mod:#9dbcf3;--mod-bg:#1b2840;--dft:#f2ad6b;--dft-bg:#3a2612;--uns:#bdb7c6;--uns-bg:#2a2730}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:var(--accent)}header,main,footer{max-width:1180px;margin:0 auto;padding:0 16px}
header nav{display:flex;flex-wrap:wrap;gap:4px 16px;padding:14px 0;border-bottom:1px solid var(--line)}
header nav a:first-child{font-weight:700;margin-right:auto;color:var(--fg);text-decoration:none}
h1{font-size:1.6rem;margin:22px 0 6px}h2{font-size:1.25rem;margin:28px 0 8px}h3{font-size:1.05rem;margin:18px 0 6px}
.lead,.muted,small{color:var(--muted)}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:.88rem;margin:6px 0 14px}
caption{text-align:left;font-weight:600;padding:4px 0;color:var(--muted)}
th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}
th{font-size:.75rem;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}td.n{white-space:nowrap;font-variant-numeric:tabular-nums}
q{font-style:italic}q::before{content:"\\201C"}q::after{content:"\\201D"}
.b{display:inline-block;font-size:.75rem;font-weight:600;padding:0 7px;border-radius:999px;white-space:nowrap}
.b-experimental{background:var(--exp-bg);color:var(--exp)}.b-model{background:var(--mod-bg);color:var(--mod)}
.b-dft{background:var(--dft-bg);color:var(--dft)}.b-unspecified{background:var(--uns-bg);color:var(--uns)}
code{font-size:.85em;overflow-wrap:anywhere}ol.src li{margin:3px 0}
.cite{border:1px solid var(--line);border-left:4px solid var(--accent);background:var(--card);border-radius:8px;padding:4px 16px 10px;margin:28px 0 0}
.cite h2{margin:10px 0 6px}.cite ul{padding-left:20px;margin:6px 0}.cite li{margin:4px 0;overflow-wrap:anywhere}
.cite pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;background:var(--row);border-radius:6px;padding:8px 10px;margin:6px 0}
footer{border-top:1px solid var(--line);margin-top:40px;padding:18px 16px 40px;color:var(--muted);font-size:.9rem}"""


# ---------------------------------------------------------------------------------------------- helpers --
def load_parity_module():
    spec = importlib.util.spec_from_file_location("build_dft_parity", os.path.join(SITE, "tools", "build_dft_parity.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def jload(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def cap(e):
    return e[:1] + e[1:].lower()


def sys_label(s):
    return "–".join(cap(e) for e in s.split("-"))


def fmt_num(v):
    if v is None or isinstance(v, bool):
        return ""
    if not isinstance(v, (int, float)):
        return str(v)
    a = abs(v)
    if a != 0 and (a < 1e-3 or a >= 1e7):
        m, e = ("%.3e" % v).split("e")
        return "%se%d" % (m, int(e))
    s = repr(float("%.6g" % v))
    return s[:-2] if s.endswith(".0") else s


def prop_name(p):
    if p in PROPN:
        return PROPN[p]
    if p and p.startswith("ACR_"):
        return "Activity " + cap(p[4:])
    return (p or "").replace("_", " ").lower()


def short_quote(q):
    q = " ".join(str(q).split())
    return q if len(q) <= QUOTE_MAX else q[:QUOTE_MAX - 1].rstrip() + "…"


def paper_short(p, pid):
    if not p:
        return pid or ""
    if p.get("authors"):
        a = p["authors"][0].split()[-1] if p["authors"][0].split() else p["authors"][0]
    else:
        a = re.split(r"[,.]", p.get("citation") or pid or "")[0].strip()
    return a + (" %s" % p["year"] if p.get("year") else "")


def doi_url(doi):
    return "https://doi.org/" + doi


def comp_text(r):
    out = []
    for pc in r.get("phase_compositions") or []:
        c = ", ".join("%s %s" % (cap(e), fmt_num(x)) for e, x in (pc.get("composition") or {}).items())
        out.append((pc["phase"] + ": " if pc.get("phase") else "") + c +
                   (" (basis unresolved)" if pc.get("basis") == "unresolved" else ""))
    if not out and r.get("overall_composition"):
        out.append("alloy: " + ", ".join("%s %s" % (cap(e), fmt_num(x)) for e, x in r["overall_composition"].items()))
    return "; ".join(out)


def phase_text(r):
    rx = r.get("reaction") or {}
    if rx.get("printed") or rx.get("type"):
        return (rx["type"] + ": " if rx.get("type") and rx.get("printed") else "") + (rx.get("printed") or rx.get("type"))
    return " + ".join(p.get("label", "") for p in r.get("phases") or [])


def value_text(r):
    if r.get("value") is None:
        if r.get("value_expression"):
            return "expression: " + short_quote(r["value_expression"])
        if r.get("value_range"):
            vr = r["value_range"]
            return "%s to %s %s" % (fmt_num(vr.get("min")), fmt_num(vr.get("max")), r.get("unit") or "")
        return ""
    s = fmt_num(r["value"]) + (" " + r["unit"] if r.get("unit") else "")
    if r.get("value_range"):
        s += " (range %s to %s)" % (fmt_num(r["value_range"].get("min")), fmt_num(r["value_range"].get("max")))
    u = r.get("uncertainty") or {}
    if u.get("basis") == "printed" and u.get("value") is not None:
        s += " ± " + fmt_num(u["value"])
    return s


def t_text(r):
    if r.get("temperature_K") is None:
        return ""
    tb = r.get("temperature_bound")
    tb = (tb.get("op") if isinstance(tb, dict) else tb) if tb else None
    return fmt_num(r["temperature_K"]) + (" (%s)" % tb if tb else "")


def jsonld_dump(obj):
    s = json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return s.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


# ------------------------------------------------------------------ document model -> HTML / Markdown --
# A cell or paragraph is a list of parts: str | ("a", text, url) | ("b", evidence_class) | ("q", quote)
# | ("code", text) | ("small", text) | ("br",). Blocks: ("p", parts) | ("h2", id, text) | ("h3", text)
# | ("ul", [parts]) | ("ol", [parts]) | ("table", caption, headers, rows, numeric_cols).
def rel(url, page_dir):
    if url.startswith(DATA_URL):
        target = url[len(DATA_URL):] or "index.html"
        r = posixpath.relpath(target, page_dir or ".")
        return r
    return url


def h_parts(parts, page_dir):
    out = []
    for p in parts:
        if isinstance(p, str):
            out.append(html.escape(p))
        elif p[0] == "a":
            out.append('<a href="%s">%s</a>' % (html.escape(rel(p[2], page_dir)), html.escape(p[1])))
        elif p[0] == "ref":
            out.append('<a href="#s%d">[%d]</a>' % (p[1], p[1]))
        elif p[0] == "b":
            out.append('<span class="b b-%s">%s</span>' % (html.escape(p[1]), html.escape(EVL.get(p[1], p[1]))))
        elif p[0] == "q":
            out.append("<q>%s</q>" % html.escape(p[1]))
        elif p[0] == "code":
            out.append("<code>%s</code>" % html.escape(p[1]))
        elif p[0] == "small":
            out.append("<small>%s</small>" % html.escape(p[1]))
        elif p[0] == "br":
            out.append("<br>")
    return "".join(out)


def md_esc(s):
    s = re.sub(r"\s+", " ", str(s))
    s = re.sub(r"([\\`*\[\]<|])", r"\\\1", s)
    return re.sub(r"(?<![A-Za-z0-9])_|_(?![A-Za-z0-9])", r"\\_", s)


def m_parts(parts, in_table=False):
    out = []
    for p in parts:
        if isinstance(p, str):
            out.append(md_esc(p))
        elif p[0] == "a":
            out.append("[%s](%s)" % (md_esc(p[1]), p[2].replace(" ", "%20").replace("(", "%28").replace(")", "%29")))
        elif p[0] == "b":
            out.append(EVL.get(p[1], p[1]))
        elif p[0] == "ref":
            out.append("[%d]" % p[1])
        elif p[0] == "q":
            out.append("“%s”" % md_esc(p[1]))
        elif p[0] == "code":
            out.append("`%s`" % p[1].replace("`", "'").replace("|", "\\|"))
        elif p[0] == "small":
            out.append("_%s_" % md_esc(p[1]))
        elif p[0] == "br":
            out.append("; " if in_table else "  \n")
    return "".join(out)


def render_html(page, page_dir):
    o = ['<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
         '<meta name="viewport" content="width=device-width,initial-scale=1">\n',
         "<title>%s</title>\n" % html.escape(page["title"]),
         '<meta name="description" content="%s">\n' % html.escape(page["description"]),
         '<link rel="canonical" href="%s">\n' % html.escape(page["url"]),
         '<link rel="alternate" type="text/markdown" href="%s">\n' % html.escape(rel(page["md_url"], page_dir)),
         '<link rel="icon" type="image/svg+xml" href="%sassets/logo/orchid-icon.svg">\n' % SITE_URL,
         "<style>%s</style>\n" % CSS]
    if page.get("jsonld"):
        o.append('<script type="application/ld+json">%s</script>\n' % jsonld_dump(page["jsonld"]))
    o.append("</head>\n<body>\n<header><nav>")
    nav = [("ORCHID data", DATA_URL + "index.html"), ("Systems", DATA_URL + "systems/index.html"),
           ("Compounds", DATA_URL + "compounds/index.html"), ("Papers", DATA_URL + "papers/index.html"),
           ("Explorer", SITE_URL + "data-explorer/"), ("orchid-project.org", SITE_URL)]
    o.append("".join('<a href="%s">%s</a>' % (html.escape(rel(u, page_dir)), html.escape(t)) for t, u in nav))
    o.append("</nav></header>\n<main>\n<h1>%s</h1>\n" % html.escape(page["h1"]))
    for b in page["blocks"]:
        k = b[0]
        if k == "p":
            o.append("<p>%s</p>\n" % h_parts(b[1], page_dir))
        elif k == "h2":
            o.append('<h2 id="%s">%s</h2>\n' % (html.escape(b[1]), html.escape(b[2])))
        elif k == "h3":
            o.append("<h3>%s</h3>\n" % html.escape(b[1]))
        elif k in ("ul", "ol"):
            o.append("<%s%s>\n" % (k, ' class="src"' if k == "ol" else ""))
            for i, it in enumerate(b[1], 1):
                o.append('<li%s>%s</li>\n' % (' id="s%d"' % i if k == "ol" else "", h_parts(it, page_dir)))
            o.append("</%s>\n" % k)
        elif k == "table":
            _, caption, headers, rows, num = b
            o.append('<div class="tw"><table>\n<caption>%s</caption>\n<thead><tr>' % html.escape(caption))
            o.append("".join('<th scope="col">%s</th>' % html.escape(t) for t in headers))
            o.append("</tr></thead>\n<tbody>\n")
            for row in rows:
                cells = []
                for i, c in enumerate(row):
                    tag = "th" if i == 0 else "td"
                    attrs = ' scope="row"' if i == 0 else (' class="n"' if i in num else "")
                    cells.append("<%s%s>%s</%s>" % (tag, attrs, h_parts(c, page_dir), tag))
                o.append("<tr>" + "".join(cells) + "</tr>\n")
            o.append("</tbody>\n</table></div>\n")
    if page.get("cite"):
        o.append(cite_html(page["cite"]))
    o.append("</main>\n<footer>%s</footer>\n</body>\n</html>\n" % h_parts(page["footer"], page_dir))
    return "".join(o)


def render_md(page):
    o = ["# %s\n\n" % md_esc(page["h1"]),
         "> %s\n\n" % md_esc(page["description"]),
         "HTML: %s · Markdown: %s\n\n" % (page["url"], page["md_url"])]
    for b in page["blocks"]:
        k = b[0]
        if k == "p":
            o.append(m_parts(b[1]) + "\n\n")
        elif k == "h2":
            o.append("## %s\n\n" % md_esc(b[2]))
        elif k == "h3":
            o.append("### %s\n\n" % md_esc(b[1]))
        elif k in ("ul", "ol"):
            for i, it in enumerate(b[1], 1):
                o.append(("%d. " % i if k == "ol" else "- ") + m_parts(it) + "\n")
            o.append("\n")
        elif k == "table":
            _, caption, headers, rows, num = b
            o.append("Table: %s\n\n" % md_esc(caption))
            o.append("| " + " | ".join(md_esc(t) for t in headers) + " |\n")
            o.append("|" + "|".join("---:" if i in num else "---" for i in range(len(headers))) + "|\n")
            for row in rows:
                o.append("| " + " | ".join(m_parts(c, True) for c in row) + " |\n")
            o.append("\n")
    if page.get("cite"):
        o.append(cite_md(page["cite"]))
    o.append("---\n\n" + m_parts(page["footer"]) + "\n")
    return "".join(o)


def cite_names(names):
    """Datasets to cite on a page, in release order, umbrella last."""
    order = [n for n in orchid_cite.NAMES if n != orchid_cite.UMBRELLA]
    return [n for n in order if n in names] + [orchid_cite.UMBRELLA]


def cite_also(names):
    seen = []
    for n in names:
        for x in orchid_cite.cite_as(n, CITE_DATA)["also_cite"]:
            if x not in seen:
                seen.append(x)
    return seen


def cite_html(names):
    names = cite_names(names)
    o = ['<aside class="cite" id="cite" aria-labelledby="cite-h">\n<h2 id="cite-h">Please cite</h2>\n<ul>\n']
    for n in names:
        e = orchid_cite.entry(n, CITE_DATA)
        o.append('<li>%s <a href="%s">doi:%s</a></li>\n' % (html.escape(e["text"].rsplit(" https://", 1)[0]),
                                                         html.escape(e["url"]), html.escape(e["doi"])))
    o.append("</ul>\n<p><small>Also cite %s. %s</small></p>\n" % (
        html.escape("; ".join(cite_also(names))), html.escape(orchid_cite.STATUS)))
    o.append("<pre>%s</pre>\n</aside>\n" % html.escape(orchid_cite.bibtex_block(names, CITE_DATA)))
    return "".join(o)


def cite_md(names):
    names = cite_names(names)
    o = ["## Please cite\n\n"]
    for n in names:
        o.append("- %s\n" % orchid_cite.entry(n, CITE_DATA)["text"])
    o.append("\nAlso cite %s.\n\n%s\n\n```bibtex\n%s\n```\n\n" % (
        "; ".join(cite_also(names)), orchid_cite.STATUS, orchid_cite.bibtex_block(names, CITE_DATA)))
    return "".join(o)


def footer_parts():
    return ["Derived data CC BY 4.0 (", ("a", "licence", LICENSE_URL), "): cite ORCHID (",
            ("a", "CITATION.cff", DATA_URL + "CITATION.cff"), ") and the original sources. Each quoted sentence is "
            "reproduced under fair use so the value can be checked; source papers, page images and full text are "
            "not redistributed. Field reference: ", ("a", "SCHEMA.md", DATA_URL + "literature/SCHEMA.md"),
            ". Machine-readable guide: ", ("a", "llms.txt", DATA_URL + "llms.txt"), "."]


def creator():
    return {"@type": "Organization", "name": "ORCHID authors", "url": SITE_URL}


def is_part_of():
    u = orchid_cite.entry(orchid_cite.UMBRELLA, CITE_DATA)
    return {"@type": "Dataset", "name": DATASET_NAME, "url": DATA_URL, "license": LICENSE_URL, "identifier": u["url"]}


# --------------------------------------------------------------------------------------------- build --
class Sources:
    """Per-page numbered source list in order of first use."""

    def __init__(self, papers, all_papers):
        self.papers, self.all, self.order = papers, all_papers, {}

    def get(self, pid):
        return self.papers.get(pid) or self.all.get(pid) or {}

    def ref(self, pid):
        if not pid:
            return []
        if pid not in self.order:
            self.order[pid] = len(self.order) + 1
        n = self.order[pid]
        return [("ref", n), " " + paper_short(self.get(pid), pid)]

    def items(self):
        out = []
        for pid, _ in sorted(self.order.items(), key=lambda kv: kv[1]):
            p = self.get(pid)
            parts = [p.get("citation") or pid]
            if p.get("doi"):
                parts += [" ", ("a", "doi:" + p["doi"], doi_url(p["doi"]))]
            out.append(parts)
        return out

    def citations(self, limit=50):
        out = []
        for pid, _ in sorted(self.order.items(), key=lambda kv: kv[1])[:limit]:
            p = self.get(pid)
            c = {"@type": "CreativeWork", "name": p.get("citation") or pid}
            if p.get("doi"):
                c["sameAs"] = doi_url(p["doi"])
                c["identifier"] = "doi:" + p["doi"]
            out.append(c)
        return out


def source_cell(r, src):
    parts = src.ref(r.get("paper"))
    origin = r.get("origin") or {}
    cs = (r.get("cited_sources") or [])
    if cs and not origin.get("own_work_of_paper"):
        c = cs[0]
        label = c.get("first_author") or re.split(r"[,.]", c.get("reference") or "")[0] or "ref. %s" % c.get("label", "?")
        label = (label + (" %s" % c["year"] if c.get("year") else "")).strip()
        parts = parts + [("br",), ("small", "cites " + label)]
        if c.get("doi"):
            parts += [" ", ("a", "doi", doi_url(c["doi"]))]
    return parts


def locator_cell(r, arts):
    art = arts.get((r.get("locator") or {}).get("ref")) or {}
    row = (r.get("locator") or {}).get("row")
    if art.get("artifact") == "table":
        s = "Table %s" % (art.get("label") or "?")
        if art.get("page") is not None:
            s += ", p. %s" % art["page"]
        if row:
            s += ", row %s" % row
        return [s]
    page = r.get("quote_page") or art.get("page")
    return ["p. %s" % page] if page is not None else []


def page_of(r, arts):
    art = arts.get((r.get("locator") or {}).get("ref")) or {}
    return r.get("quote_page") or art.get("page")


def evidence_cell(r):
    return [("b", r.get("evidence_class") or "unspecified")]


def quote_cell(r):
    return [("q", short_quote(r["quote"]))] if r.get("quote") else []


def build(data, date):
    P = load_parity_module()
    lit = os.path.join(data, "literature")
    lindex = jload(os.path.join(lit, "index.json"))
    papers_all = jload(os.path.join(lit, "papers.json"))
    parity = jload(os.path.join(lit, "dft_parity.json"))
    pindex = jload(os.path.join(data, "phase-diagrams", "index.json"))
    release = lindex.get("release", "")
    out = {}  # relative path -> text

    # ---- DFT parity, by literature record
    dft_by_rid = defaultdict(dict)
    for pt in parity["points"]:
        dft_by_rid[pt["rid"]][pt["src"]] = pt
    par_rec = parity["records"]

    # ---- calculated phase diagrams
    pd_by_sys = pindex["systems"]
    pd_docs = {}
    for s, lst in pd_by_sys.items():
        for d in lst:
            fp = os.path.join(data, "phase-diagrams", "sgte", d["file"])
            if os.path.exists(fp):
                pd_docs[d["file"]] = jload(fp)

    def hull_links(s):
        els = s.split("-")
        if len(els) != 2:
            return []
        out_ = []
        for layer, label in (("dft-mp", "Materials Project DFT hull"), ("dft-oqmd", "OQMD DFT hull")):
            rp = "hulls/%s/%s/%s.json" % (layer, els[0], s)
            if os.path.exists(os.path.join(data, rp)):
                out_.append((label, DATA_URL + rp))
        cdir = os.path.join(data, "hulls", "calphad-sgte", els[0])
        if os.path.isdir(cdir):
            for f in sorted(os.listdir(cdir)):
                if f.startswith(s + "--") and f.endswith(".json"):
                    out_.append(("CALPHAD hull " + f[len(s) + 2:-5], DATA_URL + "hulls/calphad-sgte/%s/%s" % (els[0], f)))
        return out_

    systems = sorted(set(lindex["systems"]) | set(pd_by_sys))
    sys_summary = {}
    paper_values = Counter()
    paper_systems = defaultdict(set)
    compound_lit = defaultdict(list)  # formula -> [(record, system, papers, artifacts)]
    capped = []

    for s in systems:
        meta = lindex["systems"].get(s)
        recs, arts, models, spapers = [], {}, {}, {}
        if meta:
            for f in meta["files"]:
                d = jload(os.path.join(lit, f))
                recs.extend(d["records"])
                arts.update(d.get("artifacts") or {})
                models.update(d.get("models") or {})
                spapers.update(d.get("papers") or {})
        for r in recs:
            paper_values[r.get("paper")] += 1
            paper_systems[r.get("paper")].add(s)
        acc = [r for r in recs if (r.get("flags") or {}).get("extraction_accepted")]
        n_not = len(recs) - len(acc)
        src = Sources(spapers, papers_all)
        url = DATA_URL + "systems/%s.html" % s
        md_url = DATA_URL + "systems/%s.md" % s
        label = sys_label(s)
        pds = [d for d in pd_by_sys.get(s, []) if d["file"] in pd_docs]
        blocks = []

        ev = Counter(r.get("evidence_class") for r in acc)
        pg = Counter(r.get("property_group") for r in acc)
        papers_n = len({r.get("paper") for r in acc})
        summ = "%s accepted literature values from %d paper(s)" % (format(len(acc), ","), papers_n)
        if acc:
            summ += " (" + ", ".join("%s %d" % (EVL[k], ev[k]) for k in EV_ORDER if ev.get(k)) + ")"
        summ += "; %d calculated phase diagram(s)." % len(pds)
        blocks.append(("p", [summ + (" %d further extracted value(s) not accepted by the reader-agreement checks are "
                                     "only in the JSON/CSV." % n_not if n_not else "")]))
        links = []
        if meta:
            for i, f in enumerate(meta["files"], 1):
                links += [("a", "JSON" + (" part %d" % i if len(meta["files"]) > 1 else ""), DATA_URL + "literature/" + f), " · "]
            links += [("a", "CSV", DATA_URL + "systems/%s.csv" % s), " · ",
                      ("a", "interactive explorer", SITE_URL + "data-explorer/#" + s), " · "]
        if pds:
            links += [("a", "interactive phase diagram", SITE_URL + "phase-diagrams/#%s/%s" % (s, pds[0]["file"])), " · "]
        links += [("a", "Markdown", md_url)]
        blocks.append(("p", ["Data: "] + links))
        if acc:
            rows = [[[EVL[k]], [str(ev[k])]] for k in EV_ORDER if ev.get(k)]
            blocks.append(("table", "Accepted values by evidence class", ["Evidence class", "Values"], rows, {1}))
            rows = [[[g], [str(pg[g])]] for g in sorted(pg, key=lambda g: (GROUP_ORDER.index(g) if g in GROUP_ORDER else 99, g))]
            blocks.append(("table", "Accepted values by property group", ["Property group", "Values"], rows, {1}))

        # calculated phase diagrams
        if pds:
            blocks.append(("h2", "calculated", "Calculated phase diagrams"))
            blocks.append(("p", ["Binary phase diagrams calculated with pycalphad 0.11.2 from the SGTE binary collection "
                                 "(Hallstedt 2025, CC BY 4.0). Evidence class: model."]))
            rows = []
            for d in pds:
                c = [d.get("citation") or ""]
                if d.get("doi"):
                    c += [" ", ("a", "doi:" + d["doi"], doi_url(d["doi"]))]
                rows.append([[d["tdb"]], c, [("code", d["sha256"])],
                             [("a", "viewer", SITE_URL + "phase-diagrams/#%s/%s" % (s, d["file"])), " · ",
                              ("a", "TDB", SITE_URL + "data/databases/sgte-binary/" + d["tdb"]), " · ",
                              ("a", "JSON", DATA_URL + "phase-diagrams/sgte/" + d["file"])]])
            blocks.append(("table", "Assessments with a calculated phase diagram", ["TDB file", "Assessment", "SHA-256", "Links"], rows, set()))
            rows = []
            for d in pds:
                doc = pd_docs[d["file"]]
                xc = cap(doc.get("x_component") or "")
                for inv in sorted(doc.get("invariants") or [], key=lambda i: -i["T"]):
                    pts = " + ".join("%s (x(%s) %s)" % (p_["phase"], xc, fmt_num(p_["x"])) for p_ in inv.get("points", []))
                    rows.append([[d["tdb"]], [fmt_num(inv["T"])], [pts]])
            if rows:
                blocks.append(("table", "Calculated three-phase invariants (model, mole fraction of the second component)",
                               ["TDB file", "T (K)", "Phases at the invariant"], rows, {1}))
        hl = hull_links(s)
        if hl:
            blocks.append(("p", ["Convex hulls: "] + sum([[("a", t, u), " · "] for t, u in hl], [])[:-1]))

        # literature sections
        big = len(acc) > CAP_VALUES
        if big:
            capped.append((s, len(acc)))
        inv = [r for r in acc if r.get("property_group") == "invariant reaction"]
        hform = [r for r in acc if r.get("property") == "HM_FORM"]
        rest = [r for r in acc if r.get("property_group") != "invariant reaction" and r.get("property") != "HM_FORM"]

        if inv:
            blocks.append(("h2", "invariants", "Invariant reactions (literature)"))
            rows = []
            for r in sorted(inv, key=lambda r: (-(r.get("temperature_K") or -1), r.get("id"))):
                rows.append([[t_text(r)], [phase_text(r)], [comp_text(r)], evidence_cell(r), source_cell(r, src),
                             locator_cell(r, arts), quote_cell(r)])
            blocks.append(("table", "Invariant reactions reported in the literature for the %s system (%d)" % (label, len(inv)),
                           ["T (K)", "Reaction / phases", "Compositions (mole fraction)", "Class", "Source", "Location", "Quote"],
                           rows, {0}))
        if hform:
            blocks.append(("h2", "formation", "Formation enthalpies"))
            blocks.append(("p", ["Literature values in kJ/mol-atom where the per-atom basis is known. DFT columns: the "
                                 "lowest-energy entry at the same composition (Materials Project 2026.04.13, OQMD v1.8; "
                                 "different energy references, never averaged), from ",
                                 ("a", "dft_parity.json", DATA_URL + "literature/dft_parity.json"), "."]))
            rows = []
            for r in sorted(hform, key=lambda r: (phase_text(r), r.get("value") if r.get("value") is not None else 0, r["id"])):
                cells = []
                for k in ("mp", "oqmd"):
                    pt = dft_by_rid.get(r["id"], {}).get(k)
                    cells.append([fmt_num(pt["dft"]) + " (%s)" % pt["mid"]] if pt else [])
                pr = par_rec.get(r["id"])
                comp = [("a", pr["formula"], DATA_URL + "compounds/%s.html" % pr["formula"])] if pr else []
                rows.append([[phase_text(r) or "—"], comp, [value_text(r)], [t_text(r)], evidence_cell(r)] + cells +
                            [source_cell(r, src), locator_cell(r, arts), quote_cell(r)])
            blocks.append(("table", "Enthalpies of formation, %s system (%d)" % (label, len(hform)),
                           ["Phase", "Compound", "Value", "T (K)", "Class", "MP DFT (kJ/mol-atom)", "OQMD DFT (kJ/mol-atom)",
                            "Source", "Location", "Quote"], rows, {2, 3, 5, 6}))
        if rest:
            blocks.append(("h2", "properties", "Other properties"))
            if big:
                rows = []
                cnt = Counter((r.get("property_group"), prop_name(r.get("property"))) for r in rest)
                for (g, p), n in sorted(cnt.items(), key=lambda kv: (GROUP_ORDER.index(kv[0][0]) if kv[0][0] in GROUP_ORDER else 99, kv[0])):
                    rows.append([[g or ""], [p], [str(n)]])
                blocks.append(("p", ["This system has %s accepted values; to keep the page small the %s values below are "
                                     "summarised. All values with quotes and sources are in the " % (format(len(acc), ","), format(len(rest), ","))]
                               + [("a", "JSON", DATA_URL + "literature/" + meta["files"][0]), ", the ",
                                  ("a", "CSV", DATA_URL + "systems/%s.csv" % s), " and the ",
                                  ("a", "interactive explorer", SITE_URL + "data-explorer/#" + s), "."]))
                blocks.append(("table", "Other accepted values by property (%d)" % len(rest), ["Group", "Property", "Values"], rows, {2}))
            else:
                bygroup = defaultdict(list)
                for r in rest:
                    bygroup[r.get("property_group") or "other property"].append(r)
                for g in sorted(bygroup, key=lambda g: (GROUP_ORDER.index(g) if g in GROUP_ORDER else 99, g)):
                    rs = sorted(bygroup[g], key=lambda r: (prop_name(r.get("property")), r.get("temperature_K") or 0, phase_text(r), r["id"]))
                    blocks.append(("h3", "%s (%d)" % (g[:1].upper() + g[1:], len(rs))))
                    rows = []
                    for r in rs:
                        pn = prop_name(r.get("property"))
                        if r.get("property_printed") and r["property_printed"] != pn:
                            pn = [pn, ("br",), ("small", "printed: " + short_quote(r["property_printed"]))]
                        else:
                            pn = [pn]
                        ph = phase_text(r)
                        ct = comp_text(r)
                        pc = [ph] + ([("br",), ("small", ct)] if ct else []) if ph else ([ct] if ct else [])
                        rows.append([pn, [value_text(r)], [t_text(r)], pc, evidence_cell(r), source_cell(r, src),
                                     locator_cell(r, arts), quote_cell(r)])
                    blocks.append(("table", "%s, %s system" % (g[:1].upper() + g[1:], label),
                                   ["Property", "Value", "T (K)", "Phases / composition", "Class", "Source", "Location", "Quote"],
                                   rows, {1, 2}))
        if src.order:
            blocks.append(("h2", "sources", "Sources"))
            blocks.append(("p", ["Reporting papers the values were read from (a value cited from another work names that work in its row)."]))
            blocks.append(("ol", src.items()))

        desc = "%s system: %s accepted literature values (%s) and %d calculated phase diagram(s) from ORCHID, each value with evidence class, source, page and quote." % (
            label, format(len(acc), ","), ", ".join("%s %d" % (g, pg[g]) for g in sorted(pg, key=lambda g: -pg[g])[:4]) or "none", len(pds))
        dist = []
        if meta:
            for f in meta["files"]:
                dist.append({"@type": "DataDownload", "encodingFormat": "application/json", "contentUrl": DATA_URL + "literature/" + f})
            dist.append({"@type": "DataDownload", "encodingFormat": "text/csv", "contentUrl": DATA_URL + "systems/%s.csv" % s})
        dist.append({"@type": "DataDownload", "encodingFormat": "text/markdown", "contentUrl": md_url})
        for d in pds:
            dist.append({"@type": "DataDownload", "encodingFormat": "application/json", "contentUrl": DATA_URL + "phase-diagrams/sgte/" + d["file"],
                         "name": "calculated phase diagram " + d["tdb"]})
        vm = sorted({prop_name(r.get("property")) for r in acc}) + (["Calculated phase diagram"] if pds else [])
        cits = src.citations()
        for d in pds:
            if len(cits) >= 50:
                break
            c = {"@type": "CreativeWork", "name": d.get("citation") or d["tdb"]}
            if d.get("doi"):
                c["sameAs"] = doi_url(d["doi"])
            cits.append(c)
        title = "%s system — ORCHID" % label
        jl = {"@context": "https://schema.org", "@type": "Dataset", "name": title, "description": desc, "url": url,
              "license": LICENSE_URL, "creator": creator(), "isPartOf": is_part_of(), "distribution": dist,
              "variableMeasured": vm, "citation": cits, "keywords": ["CALPHAD", "phase diagram", "thermodynamics", label],
              "isAccessibleForFree": True, "version": release, "dateModified": date}
        page = {"title": title, "h1": "%s system" % label, "description": desc, "url": url, "md_url": md_url,
                "jsonld": jl, "blocks": blocks, "footer": footer_parts(),
                "cite": (["literature"] if recs else []) + (["phase-diagrams"] if pds else []) +
                        (["hulls"] if hull_links(s) else [])}
        out["systems/%s.html" % s] = render_html(page, "systems")
        out["systems/%s.md" % s] = render_md(page)
        if meta:
            out["systems/%s.csv" % s] = system_csv(s, recs, arts, spapers, papers_all)
        sys_summary[s] = {"n": len(acc), "ev": ev, "pg": pg, "papers": papers_n, "pds": len(pds),
                          "inv": len(inv), "hform": len(hform)}

        # compounds: every HM_FORM record (any class), formula as in dft_parity
        for r in recs:
            if r.get("property") != "HM_FORM" or not (r.get("flags") or {}).get("extraction_accepted"):
                continue
            f = (par_rec.get(r["id"]) or {}).get("formula") or formula_of(P, r)
            if f:
                compound_lit[f].append((r, s, spapers, arts))

    # ---- compounds
    dft_by_formula = defaultdict(dict)
    for pt in parity["points"]:
        f = par_rec[pt["rid"]]["formula"]
        dft_by_formula[f].setdefault((pt["src"], pt["mid"]), pt)
    comp_summary = {}
    seen_lower = {}
    for f in sorted(compound_lit):
        if f.lower() in seen_lower:
            raise SystemExit("compound file name collision (case-insensitive): %s / %s" % (f, seen_lower[f.lower()]))
        seen_lower[f.lower()] = f
        items = compound_lit[f]
        url = DATA_URL + "compounds/%s.html" % f
        md_url = DATA_URL + "compounds/%s.md" % f
        src = Sources({}, papers_all)
        for _, _, sp, _ in items:
            src.papers.update(sp)
        blocks = []
        ev = Counter(r.get("evidence_class") for r, *_ in items)
        syss = sorted({s for _, s, _, _ in items})
        kj = [r["value"] for r, *_ in items if r.get("unit") == "kJ/mol-atom" and isinstance(r.get("value"), (int, float))
              and r.get("evidence_class") == "experimental"]
        dfts = dft_by_formula.get(f, {})
        summ = "%d accepted literature formation-enthalpy value(s) (%s)" % (
            len(items), ", ".join("%s %d" % (EVL[k], ev[k]) for k in EV_ORDER if ev.get(k)))
        if kj:
            kj_s = sorted(kj)
            med = kj_s[len(kj_s) // 2] if len(kj_s) % 2 else (kj_s[len(kj_s) // 2 - 1] + kj_s[len(kj_s) // 2]) / 2
            summ += "; experimental median %s kJ/mol-atom (n = %d)" % (fmt_num(round(med, 3)), len(kj))
        for k, lab in (("mp", "Materials Project"), ("oqmd", "OQMD")):
            vals = [pt for (sk, _), pt in dfts.items() if sk == k]
            if vals:
                best = min(vals, key=lambda p_: (p_["dft"], p_["mid"]))
                summ += "; %s %s kJ/mol-atom (%s)" % (lab, fmt_num(best["dft"]), best["mid"])
        blocks.append(("p", [summ + "."]))
        blocks.append(("p", ["Systems: "] + sum([[("a", sys_label(s), DATA_URL + "systems/%s.html" % s), " · "] for s in syss], [])[:-1]))
        rows = []
        for r, s, sp, arts in sorted(items, key=lambda it: (EV_ORDER.index(it[0].get("evidence_class")) if it[0].get("evidence_class") in EV_ORDER else 9,
                                                               it[0].get("value") if isinstance(it[0].get("value"), (int, float)) else 0, it[0]["id"])):
            rows.append([[value_text(r)], [phase_text(r)], [t_text(r)], evidence_cell(r),
                         [r.get("method") or ""], [("a", sys_label(s), DATA_URL + "systems/%s.html" % s)],
                         source_cell(r, src), locator_cell(r, arts), quote_cell(r)])
        blocks.append(("h2", "literature", "Literature values"))
        blocks.append(("table", "Enthalpy of formation of %s reported in the literature (%d)" % (f, len(items)),
                       ["Value", "Phase as printed", "T (K)", "Class", "Method", "System", "Source", "Location", "Quote"],
                       rows, {0, 2}))
        if dfts:
            blocks.append(("h2", "dft", "DFT formation energies"))
            blocks.append(("p", ["Lowest-energy entries at this composition. Materials Project and OQMD use different "
                                 "energy references and are never averaged. 1 eV/atom = 96.485 kJ/mol-atom."]))
            rows = []
            for (k, mid), pt in sorted(dfts.items()):
                rows.append([[parity["dft_sources"][k]["label"]], [mid], [fmt_num(pt["dft"])], [fmt_num(pt["dft_eV"])],
                             [pt.get("sg") or ""], [fmt_num(pt.get("ehull"))], ["yes" if pt.get("theoretical") else "no"]])
            blocks.append(("table", "DFT formation energy of %s" % f,
                           ["Source", "Entry", "kJ/mol-atom", "eV/atom", "Space group", "E above hull (eV/atom)", "Theoretical"],
                           rows, {2, 3, 5}))
            blocks.append(("p", ["DFT sources: "] + sum([[parity["dft_sources"][k]["citation"], " ",
                                                          ("a", "doi:" + parity["dft_sources"][k]["doi"], doi_url(parity["dft_sources"][k]["doi"])),
                                                          " (%s); " % parity["dft_sources"][k]["license"]] for k in ("mp", "oqmd")], [])))
        blocks.append(("h2", "sources", "Sources"))
        blocks.append(("ol", src.items()))
        desc = "Enthalpy of formation of %s: %d literature value(s) with sources, pages and quotes%s, from ORCHID." % (
            f, len(items), " and DFT values (Materials Project / OQMD)" if dfts else "")
        jl = {"@context": "https://schema.org", "@type": "Dataset", "name": "%s formation enthalpy — ORCHID" % f,
              "description": desc, "url": url, "license": LICENSE_URL, "creator": creator(), "isPartOf": is_part_of(),
              "distribution": [{"@type": "DataDownload", "encodingFormat": "text/markdown", "contentUrl": md_url},
                               {"@type": "DataDownload", "encodingFormat": "application/json",
                                "contentUrl": DATA_URL + "literature/dft_parity.json"}],
              "variableMeasured": ["Enthalpy of formation"] + (["DFT formation energy"] if dfts else []),
              "citation": src.citations(), "keywords": ["formation enthalpy", f, "CALPHAD", "DFT"],
              "isAccessibleForFree": True, "version": release, "dateModified": date}
        page = {"title": "%s formation enthalpy — ORCHID" % f, "h1": "%s: enthalpy of formation" % f, "description": desc,
                "url": url, "md_url": md_url, "jsonld": jl, "blocks": blocks, "footer": footer_parts(),
                "cite": ["literature", "dft-benchmark"]}
        out["compounds/%s.html" % f] = render_html(page, "compounds")
        out["compounds/%s.md" % f] = render_md(page)
        comp_summary[f] = {"n": len(items), "ev": ev, "systems": syss,
                           "mp": min((pt["dft"] for (k, _), pt in dfts.items() if k == "mp"), default=None),
                           "oqmd": min((pt["dft"] for (k, _), pt in dfts.items() if k == "oqmd"), default=None),
                           "exp_n": len(kj)}

    # ---- indexes
    index_pages(out, sys_summary, comp_summary, papers_all, paper_values, paper_systems, lindex, pindex, date, release)
    return out, capped


def formula_of(P, r):
    system = [e.capitalize() for e in (r.get("elements") or r.get("system", "").split("-")) if e]
    found = {}
    for p in r.get("phases") or []:
        if p.get("label"):
            c, _ = P.parse_phase(p["label"], system)
            if c:
                found[P.reduced_key(c)] = c
    if len(found) == 1:
        return P.formula_str(next(iter(found.values())))
    if found:
        return None
    comps = {tuple(sorted((k.capitalize(), round(float(x), 4)) for k, x in pc["composition"].items() if float(x) > 0))
             for pc in (r.get("phase_compositions") or [])
             if pc.get("basis") == "mole_fraction" and pc.get("composition_source") in ("stoichiometry", "formula_nominal")}
    if len(comps) == 1:
        key = next(iter(comps))
        if len(key) >= 2 and abs(sum(x for _, x in key) - 1) < 2e-3:
            return P.formula_str({el: Fraction(x).limit_denominator(60) for el, x in key})
    return None


def system_csv(s, recs, arts, spapers, papers_all):
    cols = ["id", "system", "property", "property_group", "property_printed", "record_kind", "value", "unit",
            "value_expression", "uncertainty", "uncertainty_basis", "temperature_K", "temperature_printed", "pressure_Pa",
            "phases", "reaction_type", "reaction_printed", "compositions", "evidence_class", "method",
            "extraction_accepted", "paper", "paper_doi", "paper_citation", "location", "page", "quote"]
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for r in recs:
        p = spapers.get(r.get("paper")) or papers_all.get(r.get("paper")) or {}
        u = r.get("uncertainty") or {}
        rx = r.get("reaction") or {}
        w.writerow([r.get("id"), s, r.get("property"), r.get("property_group"), r.get("property_printed") or "",
                    r.get("record_kind"), "" if r.get("value") is None else r["value"], r.get("unit") or "",
                    r.get("value_expression") or "", "" if u.get("value") is None else u["value"], u.get("basis") or "",
                    "" if r.get("temperature_K") is None else r["temperature_K"], r.get("temperature_printed") or "",
                    "" if r.get("pressure_Pa") is None else r["pressure_Pa"],
                    "; ".join(x.get("label", "") for x in r.get("phases") or []), rx.get("type") or "", rx.get("printed") or "",
                    comp_text(r), r.get("evidence_class"), r.get("method") or "",
                    "true" if (r.get("flags") or {}).get("extraction_accepted") else "false", r.get("paper") or "",
                    p.get("doi") or "", p.get("citation") or "", "".join(x for x in locator_cell(r, arts) if isinstance(x, str)),
                    page_of(r, arts) or "", r.get("quote") or ""])
    return buf.getvalue()


def index_pages(out, sys_summary, comp_summary, papers_all, paper_values, paper_systems, lindex, pindex, date, release):
    foot = footer_parts()
    # systems index
    rows = []
    for s in sorted(sys_summary, key=lambda k: (len(k.split("-")), k)):
        v = sys_summary[s]
        rows.append([[("a", sys_label(s), DATA_URL + "systems/%s.html" % s)], [str(v["n"])],
                     [", ".join("%s %d" % (EVL[k], v["ev"][k]) for k in EV_ORDER if v["ev"].get(k))],
                     [str(v["inv"])], [str(v["hform"])], [str(v["pds"])],
                     [("a", "md", DATA_URL + "systems/%s.md" % s)]])
    url, md_url = DATA_URL + "systems/index.html", DATA_URL + "systems/index.md"
    desc = "%d chemical systems with literature thermodynamic data and/or calculated phase diagrams (ORCHID)." % len(sys_summary)
    blocks = [("p", ["One page per chemical system (elements sorted alphabetically, joined by '-'), each with a Markdown "
                     "twin at systems/<SYS>.md and a CSV at systems/<SYS>.csv. Binary systems first."]),
              ("table", "Chemical systems (%d)" % len(rows), ["System", "Accepted values", "By evidence class", "Invariants",
                                                              "Formation enthalpies", "Calculated diagrams", "Markdown"], rows, {1, 3, 4, 5})]
    jl = {"@context": "https://schema.org", "@type": "Dataset", "name": "ORCHID chemical systems", "description": desc, "url": url,
          "license": LICENSE_URL, "creator": creator(), "isPartOf": is_part_of(), "dateModified": date,
          "hasPart": [DATA_URL + "systems/%s.html" % s for s in sorted(sys_summary)][:1000]}
    page = {"title": "Chemical systems — ORCHID", "h1": "Chemical systems", "description": desc, "url": url, "md_url": md_url,
            "jsonld": jl, "blocks": blocks, "footer": foot, "cite": ["literature", "phase-diagrams", "hulls"]}
    out["systems/index.html"] = render_html(page, "systems")
    out["systems/index.md"] = render_md(page)

    # compounds index
    rows = []
    for f in sorted(comp_summary, key=lambda k: k.lower()):
        v = comp_summary[f]
        rows.append([[("a", f, DATA_URL + "compounds/%s.html" % f)], [str(v["n"])], [str(v["exp_n"])],
                     [fmt_num(v["mp"]) if v["mp"] is not None else ""], [fmt_num(v["oqmd"]) if v["oqmd"] is not None else ""],
                     sum([[("a", sys_label(s), DATA_URL + "systems/%s.html" % s), " "] for s in v["systems"][:6]], [])])
    url, md_url = DATA_URL + "compounds/index.html", DATA_URL + "compounds/index.md"
    desc = "%d compounds with literature formation-enthalpy values, compared with Materials Project and OQMD DFT (ORCHID)." % len(comp_summary)
    blocks = [("p", ["Compound ids are reduced formulas (metals first, then alphabetical). Experimental values with the "
                     "per-atom basis known are in kJ/mol-atom; DFT columns give the lowest-energy entry at the composition."]),
              ("table", "Compounds (%d)" % len(rows), ["Compound", "Literature values", "Experimental kJ/mol-atom values",
                                                      "MP DFT (kJ/mol-atom)", "OQMD DFT (kJ/mol-atom)", "Systems"], rows, {1, 2, 3, 4})]
    jl = {"@context": "https://schema.org", "@type": "Dataset", "name": "ORCHID compound formation enthalpies", "description": desc,
          "url": url, "license": LICENSE_URL, "creator": creator(), "isPartOf": is_part_of(), "dateModified": date,
          "variableMeasured": ["Enthalpy of formation", "DFT formation energy"]}
    page = {"title": "Compounds — ORCHID", "h1": "Compounds: enthalpies of formation", "description": desc, "url": url,
            "md_url": md_url, "jsonld": jl, "blocks": blocks, "footer": foot, "cite": ["literature", "dft-benchmark"]}
    out["compounds/index.html"] = render_html(page, "compounds")
    out["compounds/index.md"] = render_md(page)

    # papers index
    rep = [pid for pid, p in papers_all.items() if "reporting_paper" in (p.get("roles") or []) and paper_values.get(pid)]
    rows = []
    for pid in sorted(rep, key=lambda k: (-(papers_all[k].get("year") or 0), (papers_all[k].get("citation") or k).lower())):
        p = papers_all[pid]
        syss = sorted(paper_systems[pid])
        rows.append([[p.get("citation") or pid], [("a", "doi:" + p["doi"], doi_url(p["doi"]))] if p.get("doi") else [pid],
                     sum([[("a", sys_label(s), DATA_URL + "systems/%s.html" % s), " "] for s in syss[:12]], []) +
                     (["… %d more" % (len(syss) - 12)] if len(syss) > 12 else []),
                     [str(paper_values[pid])]])
    url, md_url = DATA_URL + "papers/index.html", DATA_URL + "papers/index.md"
    desc = "%d published papers from which ORCHID literature values were read, with DOI, systems and number of values." % len(rows)
    blocks = [("p", ["Each listed paper is a reporting paper: values were read from its tables or running text. The papers "
                     "themselves are not redistributed; follow the DOI. Number of values counts all extracted records, "
                     "accepted or not. Full bibliography including cited works: ",
                     ("a", "papers.json", DATA_URL + "literature/papers.json"), "."]),
              ("table", "Reporting papers (%d), newest first" % len(rows), ["Citation", "DOI", "Systems", "Values"], rows, {3})]
    jl = {"@context": "https://schema.org", "@type": "Dataset", "name": "ORCHID reporting papers", "description": desc, "url": url,
          "license": LICENSE_URL, "creator": creator(), "isPartOf": is_part_of(), "dateModified": date,
          "citation": [{"@type": "CreativeWork", "name": papers_all[pid].get("citation") or pid,
                        **({"sameAs": doi_url(papers_all[pid]["doi"])} if papers_all[pid].get("doi") else {})} for pid in sorted(rep)[:50]]}
    page = {"title": "Reporting papers — ORCHID", "h1": "Reporting papers", "description": desc, "url": url, "md_url": md_url,
            "jsonld": jl, "blocks": blocks, "footer": foot, "cite": ["literature"]}
    out["papers/index.html"] = render_html(page, "papers")
    out["papers/index.md"] = render_md(page)

    # root index
    t = lindex["totals"]
    n_pd = sum(len(v) for v in pindex["systems"].values())
    url, md_url = DATA_URL, DATA_URL + "index.md"
    desc = ("ORCHID data: %s literature values (invariant reactions, phase transitions, phase compositions, thermochemical "
            "and structural data) from %s papers, %s calculated binary phase diagrams and DFT/CALPHAD convex hulls; CC BY 4.0."
            % (format(t["records"], ","), format(t["papers_reporting"], ","), format(n_pd, ",")))
    blocks = [
        ("p", ["This is the data site of ORCHID (Open Repository for CALPHAD Hypotheses, Intelligence, and Discovery, ",
               ("a", "orchid-project.org", SITE_URL), "). Every value carries its evidence class (experimental, model, "
               "DFT, unspecified; never mixed or averaged), its source paper and DOI, and the table row or one verbatim "
               "sentence and page it was read from. Release: ", ("code", release), "."]),
        ("p", ["Status: " + orchid_cite.STATUS]),
        ("h2", "browse", "Browse without JavaScript"),
        ("ul", [[("a", "Chemical systems", DATA_URL + "systems/index.html"), " — %d pages (also " % len(sys_summary),
                 ("a", "Markdown", DATA_URL + "systems/index.md"), "): literature values, invariants, formation enthalpies, calculated phase diagrams."],
                [("a", "Compounds", DATA_URL + "compounds/index.html"), " — %d formation-enthalpy pages with DFT comparison." % len(comp_summary)],
                [("a", "Reporting papers", DATA_URL + "papers/index.html"), " — citation, DOI, systems, number of values."],
                ["Interactive tools: ", ("a", "data explorer", SITE_URL + "data-explorer/"), ", ",
                 ("a", "phase diagrams", SITE_URL + "phase-diagrams/"), ", ", ("a", "DFT vs experiment", SITE_URL + "dft-parity/"), "."]]),
        ("h2", "folders", "Folder guide"),
        ("table", "Folders of the data site", ["Path", "Content"], [
            [[("a", "literature/", DATA_URL + "literature/SCHEMA.md")], ["Literature records: index.json, systems/<SYS>.json, papers.json, contexts.json, dft_parity.json, full CSV/JSONL dumps; fields in SCHEMA.md."]],
            [[("a", "phase-diagrams/", DATA_URL + "phase-diagrams/index.json")], ["Calculated binary phase diagrams (pycalphad, SGTE binary collection): index.json and sgte/<TDB>.json (phase regions, invariants, grid)."]],
            [["hulls/"], ["Binary convex hulls: dft-mp/ (Materials Project), dft-oqmd/ (OQMD), calphad-sgte/ (assessed CALPHAD), one JSON per system."]],
            [[("a", "systems/", DATA_URL + "systems/index.html")], ["Static pages per chemical system: <SYS>.html, <SYS>.md, <SYS>.csv."]],
            [[("a", "compounds/", DATA_URL + "compounds/index.html")], ["Static pages per compound (reduced formula): <FORMULA>.html, <FORMULA>.md."]],
            [[("a", "papers/", DATA_URL + "papers/index.html")], ["Index of reporting papers."]],
            [["parquet/"], ["Columnar copies of the tables, where present."]],
            [["zenodo/"], ["Zenodo release metadata: one deposition record per citable dataset plus the umbrella record, with file manifests (sizes, SHA-256)."]],
            [[("a", "llms.txt", DATA_URL + "llms.txt"), ", ", ("a", "llms-full.txt", DATA_URL + "llms-full.txt")], ["Guides for language models and agents."]],
        ], set()),
        ("h2", "licence", "Licence and citation"),
        ("p", ["Release %s is archived on Zenodo as four citable datasets, each with its own DOI, linked by an umbrella "
               "record: " % orchid_cite.RELEASE_TAG] +
              sum([[("a", orchid_cite.entry(n, CITE_DATA)["title"], orchid_cite.entry(n, CITE_DATA)["url"]),
                    "; " if i < len(orchid_cite.NAMES) - 1 else "."] for i, n in enumerate(orchid_cite.NAMES)], [])),
        ("p", ["Derived data: CC BY 4.0. Cite the ORCHID release (", ("a", "CITATION.cff", DATA_URL + "CITATION.cff"),
               ") and the original source of every value you use (its paper and DOI). Quoted sentences are single sentences "
               "reproduced under fair use so values can be checked; source papers, page images and full text are not "
               "redistributed. Licences of the DFT and CALPHAD sources: ", ("a", "LICENSES.md", DATA_URL + "LICENSES.md"),
               ", ", ("a", "NOTICE.md", DATA_URL + "NOTICE.md"), "."]),
    ]
    jl = {"@context": "https://schema.org", "@type": "Dataset", "name": DATASET_NAME, "description": desc, "url": url,
          "license": LICENSE_URL, "creator": creator(), "isAccessibleForFree": True, "version": release, "dateModified": date,
          "keywords": ["CALPHAD", "phase diagrams", "thermodynamics", "formation enthalpy", "DFT", "convex hull"],
          "distribution": [{"@type": "DataDownload", "encodingFormat": "application/gzip", "contentUrl": DATA_URL + "literature/orchid-literature.csv.gz"},
                           {"@type": "DataDownload", "encodingFormat": "application/gzip", "contentUrl": DATA_URL + "literature/orchid-literature.jsonl.gz"},
                           {"@type": "DataDownload", "encodingFormat": "application/json", "contentUrl": DATA_URL + "literature/index.json"}]}
    page = {"title": "ORCHID data", "h1": "ORCHID data", "description": desc, "url": url, "md_url": md_url, "jsonld": jl,
            "blocks": blocks, "footer": foot, "cite": list(orchid_cite.NAMES)}
    out["index.html"] = render_html(page, "")
    out["index.md"] = render_md(page)

    # llms.txt / llms-full.txt
    out["llms.txt"] = "\n".join([
        "# ORCHID data", "",
        "> %s" % desc, "",
        "Evidence classes are never mixed: experimental | model (CALPHAD, linked to the TDB by SHA-256) | dft | unspecified.",
        "Formation enthalpies are in kJ/mol-atom when the per-atom basis is known; temperatures in K; compositions in mole fraction.", "",
        "## Per system (SYS = element symbols, uppercase, sorted alphabetically, joined by '-', e.g. AL-NI)",
        "- Markdown page: %ssystems/<SYS>.md (summary, calculated diagrams, invariants, formation enthalpies with DFT, other properties, sources)" % DATA_URL,
        "- HTML page: %ssystems/<SYS>.html" % DATA_URL,
        "- CSV: %ssystems/<SYS>.csv (all records incl. not accepted; column extraction_accepted)" % DATA_URL,
        "- Full JSON records: %sliterature/systems/<SYS>.json (large systems: also <SYS>.part2.json ...; see literature/index.json)" % DATA_URL,
        "- Calculated phase diagrams: %sphase-diagrams/index.json -> phase-diagrams/sgte/<file>.json" % DATA_URL, "",
        "## Indexes",
        "- [Systems](%ssystems/index.md)" % DATA_URL,
        "- [Compounds (formation enthalpies, literature vs Materials Project / OQMD)](%scompounds/index.md); per compound: %scompounds/<FORMULA>.md" % (DATA_URL, DATA_URL),
        "- [Reporting papers](%spapers/index.md)" % DATA_URL,
        "- [All systems with one-line summaries](%sllms-full.txt)" % DATA_URL,
        "- [Field reference](%sliterature/SCHEMA.md)" % DATA_URL, "",
        "## Bulk files",
        "- %sliterature/orchid-literature.jsonl.gz (all records, one JSON object per line)" % DATA_URL,
        "- %sliterature/orchid-literature.csv.gz (all records, flattened)" % DATA_URL,
        "- %sliterature/dft_parity.json (literature formation enthalpies vs DFT)" % DATA_URL, "",
        "## Licence and citation",
        "- Status: %s" % orchid_cite.STATUS,
        "- Zenodo release %s (cite the dataset you use and the umbrella record):" % orchid_cite.RELEASE_TAG] + [
        "  - %s: https://doi.org/%s" % (orchid_cite.entry(n, CITE_DATA)["title"], orchid_cite.entry(n, CITE_DATA)["doi"])
        for n in orchid_cite.NAMES] + [
        "- Every JSON file carries a top-level cite_as object (text, BibTeX, DOI); Parquet files carry orchid.cite_as metadata; literature/CITATION.txt sits beside the CSV.",
        "- Derived data CC BY 4.0. Cite the ORCHID release (%sCITATION.cff) and each value's source paper (DOI given per value)." % DATA_URL,
        "- Quotes: one verbatim sentence per value, fair use. Source papers, page images and full text are not redistributed.",
        "- Use values with extraction_accepted = true unless you handle disputes; uncertainty with basis stated_default is a placeholder, not a measurement.",
        "- Project site: %s (interactive explorer %sdata-explorer/#<SYS>)" % (SITE_URL, SITE_URL), ""])
    lines = ["# ORCHID data: all systems", "", "> One line per chemical system: accepted literature values by evidence class, "
             "invariants, formation enthalpies, calculated phase diagrams. Append .md to read a system.", "",
             "Release: %s. Licence: CC BY 4.0 (cite ORCHID and each source paper)." % release, "", "## Systems", ""]
    for s in sorted(sys_summary, key=lambda k: (len(k.split("-")), k)):
        v = sys_summary[s]
        parts = ["%d values" % v["n"]]
        evs = ", ".join("%s %d" % (EVL[k], v["ev"][k]) for k in EV_ORDER if v["ev"].get(k))
        if evs:
            parts[0] += " (%s)" % evs
        top = ", ".join("%s %d" % (g, n) for g, n in sorted(v["pg"].items(), key=lambda kv: (-kv[1], kv[0]))[:3])
        if top:
            parts.append(top)
        if v["pds"]:
            parts.append("%d calculated diagram(s)" % v["pds"])
        lines.append("- [%s](%ssystems/%s.md): %s" % (sys_label(s), DATA_URL, s, "; ".join(parts)))
    lines += ["", "## Compounds (formation enthalpy)", ""]
    for f in sorted(comp_summary, key=lambda k: k.lower()):
        v = comp_summary[f]
        extra = []
        if v["mp"] is not None:
            extra.append("MP %s" % fmt_num(v["mp"]))
        if v["oqmd"] is not None:
            extra.append("OQMD %s" % fmt_num(v["oqmd"]))
        lines.append("- [%s](%scompounds/%s.md): %d literature value(s)%s" % (f, DATA_URL, f, v["n"],
                     ("; DFT kJ/mol-atom " + ", ".join(extra)) if extra else ""))
    out["llms-full.txt"] = "\n".join(lines) + "\n"


def robots(sitemaps):
    lines = ["# ORCHID: all crawlers, including AI crawlers, are welcome.", ""]
    for b in AI_BOTS:
        lines += ["User-agent: %s" % b, "Allow: /", ""]
    lines += ["User-agent: *", "Allow: /", ""] + ["Sitemap: %s" % u for u in sitemaps]
    return "\n".join(lines) + "\n"


def urlset(urls, date):
    o = ['<?xml version="1.0" encoding="UTF-8"?>', '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    for u in urls:
        o.append("<url><loc>%s</loc><lastmod>%s</lastmod></url>" % (html.escape(u), date))
    o.append("</urlset>")
    return "\n".join(o) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=DEFAULT_DATA, help="orchid-data repository (inputs read, pages written)")
    ap.add_argument("--site", default=SITE, help="website repository (sitemap.xml, sitemap-pages.xml, robots.txt)")
    ap.add_argument("--date", default="2026-10-22", help="lastmod date for the sitemaps (YYYY-MM-DD)")
    ap.add_argument("--check", action="store_true", help="exit 1 if outputs are missing, stale or orphaned")
    a = ap.parse_args(argv)
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", a.date):
        ap.error("--date must be YYYY-MM-DD")
    global CITE_DATA
    CITE_DATA = a.data
    data_out, capped = build(a.data, a.date)
    pages = sorted(p for p in data_out if p.endswith(".html"))
    urls = [DATA_URL + ("" if p == "index.html" else p) for p in pages]
    data_out["sitemap.xml"] = urlset(urls, a.date)
    data_out["robots.txt"] = robots([DATA_URL + "sitemap.xml"])
    site_out = {
        "sitemap-pages.xml": urlset([SITE_URL + p for p in SITE_PAGES], a.date),
        "sitemap.xml": "\n".join(['<?xml version="1.0" encoding="UTF-8"?>',
                                  '<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">',
                                  "<sitemap><loc>%ssitemap-pages.xml</loc><lastmod>%s</lastmod></sitemap>" % (SITE_URL, a.date),
                                  "<sitemap><loc>%ssitemap.xml</loc><lastmod>%s</lastmod></sitemap>" % (DATA_URL, a.date),
                                  "</sitemapindex>"]) + "\n",
        "robots.txt": robots([SITE_URL + "sitemap.xml", DATA_URL + "sitemap.xml"]),
    }
    targets = [(a.data, data_out), (a.site, site_out)]
    bad = []
    for root, files in targets:
        for rp, text in sorted(files.items()):
            fp = os.path.join(root, rp)
            b = text.encode("utf-8")
            if a.check:
                try:
                    with open(fp, "rb") as fh:
                        if fh.read() != b:
                            bad.append("stale: " + fp)
                except FileNotFoundError:
                    bad.append("missing: " + fp)
            else:
                os.makedirs(os.path.dirname(fp), exist_ok=True)
                try:
                    with open(fp, "rb") as fh:
                        same = fh.read() == b
                except FileNotFoundError:
                    same = False
                if not same:
                    with open(fp, "wb") as fh:
                        fh.write(b)
    for d in OWN_DIRS:  # orphans in our own folders
        dp = os.path.join(a.data, d)
        if os.path.isdir(dp):
            for f in sorted(os.listdir(dp)):
                rp = d + "/" + f
                if rp not in data_out:
                    if a.check:
                        bad.append("orphan: " + os.path.join(dp, f))
                    else:
                        os.remove(os.path.join(dp, f))
    n_html = len(pages)
    if a.check:
        if bad:
            print("\n".join(bad[:50]) + ("\n... %d more" % (len(bad) - 50) if len(bad) > 50 else ""))
            print("FAIL: %d problem(s)" % len(bad))
            return 1
        print("OK: %d files current (%d HTML pages)" % (sum(len(f) for _, f in targets), n_html))
        return 0
    kinds = Counter(p.split("/")[0] if "/" in p else "root" for p in data_out)
    print("wrote %d data files (%s), %d site files; %d HTML pages in the data sitemap" % (
        len(data_out), ", ".join("%s %d" % kv for kv in sorted(kinds.items())), len(site_out), len(urls)))
    for s, n in capped:
        print("capped: %s (%d accepted values; non-invariant, non-HM_FORM values summarised)" % (s, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
