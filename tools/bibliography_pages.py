#!/usr/bin/env python3
"""Static, JavaScript-free HTML pages of the ORCHID bibliography (hook for tools/build_static_pages.py).

  bibliography/pages/papers/<key>.html    one page per paper: metadata, roles, evidence counts, relations
  bibliography/pages/systems/<SYS>.html   one page per chemical system: assessments and the data papers they used
  bibliography/pages/index.html           entry page (counts, links to system pages)

Every paper page embeds schema.org JSON-LD (ScholarlyArticle) and a COinS span (Z3988) so reference managers
(Zotero, Mendeley) and crawlers recognise the work described. Highwire `citation_*` tags are deliberately not used:
they tell Google Scholar that the page *is* the article (hosting its full text), which these pages are not.

Called by tools/build_bibliography.py (`render`, pages included in its --check) and usable as a hook by
tools/build_static_pages.py:

  import bibliography_pages
  out.update(bibliography_pages.pages_from_files(DATA))   # {"bibliography/pages/...": html}

Standard library only.
"""
from __future__ import annotations

import html
import json
import os
import urllib.parse

DATA_URL = "https://data.orchid-project.org/"
REL_MAX = 50  # relations listed per group on a page (all of them are in papers/<key>.json)
STYLE = """:root{--bg:#fbfaf7;--fg:#1d1b20;--muted:#5e5a66;--line:#e2dfe6;--card:#fff;--accent:#7b3fa0;--row:#f6f4f8;
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
code{font-size:.85em;overflow-wrap:anywhere}main p,main li{overflow-wrap:break-word}
pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;background:var(--row);border-radius:6px;padding:8px 10px;margin:6px 0}ol.src li{margin:3px 0}
.cite{border:1px solid var(--line);border-left:4px solid var(--accent);background:var(--card);border-radius:8px;padding:4px 16px 10px;margin:28px 0 0}
.cite h2{margin:10px 0 6px}.cite ul{padding-left:20px;margin:6px 0}.cite li{margin:4px 0;overflow-wrap:anywhere}
.cite pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;background:var(--row);border-radius:6px;padding:8px 10px;margin:6px 0}
footer{border-top:1px solid var(--line);margin-top:40px;padding:18px 16px 40px;color:var(--muted);font-size:.9rem}"""
ROLE_LABEL = {
    "corpus_assessment_read": "read by ORCHID (values extracted)",
    "cited_original_source": "cited as the original source of values",
    "tdb_source_assessment": "assessment behind a TDB file",
    "tdb_compilation": "TDB collection (compilation)",
    "handbook_diagram_source": "original evaluation behind an assessed (handbook) diagram",
    "discovery_lead": "discovery lead (found by Crossref assessment discovery, not read)",
    "context_source": "same-system context source",
    "espei_reference": "cited by an ESPEI dataset",
}
EDGE_LABEL = {
    ("cites", "out"): "Cites (values carry this citation)",
    ("cites", "in"): "Cited for values by",
    ("used_in_assessment", "out"): "Used in (assessments / papers that took values from it)",
    ("used_in_assessment", "in"): "Data papers used",
    ("models_system_in_tdb", "out"): "TDB files modelling its system",
    ("reassesses", "out"): "Reassesses (cites an earlier assessment of the same system)",
    ("reassesses", "in"): "Reassessed by",
    ("handbook_reproduces", "in"): "Reproduced by handbook diagrams",
    ("cited_by_crossref", "out"): "References (Crossref, publisher-deposited)",
    ("cited_by_crossref", "in"): "Cited by (Crossref reference lists)",
}


def esc(s):
    return html.escape("" if s is None else str(s), quote=True)


def sys_label(s):
    return "–".join(e.capitalize() for e in s.split("-"))


def author_text(a):
    if a.get("family"):
        return (a.get("given") + " " if a.get("given") else "") + a["family"]
    return a.get("name") or ""


def coins(r):
    """OpenURL 1.0 ContextObject in a span (COinS)."""
    kv = [("ctx_ver", "Z39.88-2004")]
    if r.get("type") in (None, "journal-article"):
        kv += [("rft_val_fmt", "info:ofi/fmt:kev:mtx:journal"), ("rft.genre", "article")]
        if r.get("title"):
            kv.append(("rft.atitle", r["title"]))
        if r.get("container_title"):
            kv.append(("rft.jtitle", r["container_title"]))
    else:
        kv += [("rft_val_fmt", "info:ofi/fmt:kev:mtx:book"), ("rft.genre", "book" if "book" in (r.get("type") or "") else "unknown")]
        if r.get("title"):
            kv.append(("rft.btitle" if r.get("type") in ("book", "monograph", "edited-book") else "rft.atitle", r["title"]))
        if r.get("container_title"):
            kv.append(("rft.btitle", r["container_title"]))
    for a in r.get("authors") or []:
        if a.get("family"):
            kv.append(("rft.au", (a["family"] + (", " + a["given"] if a.get("given") else ""))))
        elif a.get("name"):
            kv.append(("rft.au", a["name"]))
    for k, f in (("year", "rft.date"), ("volume", "rft.volume"), ("issue", "rft.issue"), ("pages", "rft.pages")):
        if r.get(k):
            kv.append((f, str(r[k])))
    for i in r.get("issn") or []:
        kv.append(("rft.issn", i))
    if r.get("doi"):
        kv.append(("rft_id", "info:doi/" + r["doi"]))
        kv.append(("rft_id", "https://doi.org/" + r["doi"]))
    return '<span class="Z3988" title="%s"></span>' % esc(urllib.parse.urlencode(kv))


def page(title, desc, canonical, body, jsonld=None, depth=2):
    up = "../" * depth
    ld = ""
    if jsonld is not None:
        ld = '<script type="application/ld+json">%s</script>\n' % json.dumps(jsonld, ensure_ascii=False, sort_keys=True,
                                                                         separators=(",", ":")).replace("</", "<\\/")
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n"
        f"<title>{esc(title)}</title>\n<meta name=\"description\" content=\"{esc(desc)}\">\n"
        f"<link rel=\"canonical\" href=\"{esc(canonical)}\">\n"
        "<link rel=\"icon\" type=\"image/svg+xml\" href=\"https://orchid-project.org/assets/logo/orchid-icon.svg\">\n"
        f"<link rel=\"stylesheet\" href=\"{'../' if depth == 2 else ''}style.css\">\n{ld}</head>\n<body>\n"
        f"<header><nav><a href=\"{up}../index.html\">ORCHID data</a><a href=\"{up}pages/index.html\">Bibliography</a>"
        f"<a href=\"{up}../systems/index.html\">Systems</a><a href=\"{up}../papers/index.html\">Papers (values)</a>"
        "<a href=\"https://orchid-project.org/\">orchid-project.org</a></nav></header>\n<main>\n"
        + body +
        "</main>\n<footer>Bibliographic metadata from Crossref (CC0); ORCHID-derived fields (roles, evidence counts, "
        "relations) CC BY 4.0. ORCHID does not redistribute the papers. Field reference: "
        f"<a href=\"{up}SCHEMA.md\">SCHEMA.md</a>.</footer>\n</body>\n</html>\n")


def _counter_table(caption, d, limit=None):
    items = sorted(d.items(), key=lambda kv: (-kv[1], kv[0]))
    items = items[:limit] if limit else items
    rows = "".join(f"<tr><td>{esc(k)}</td><td class=\"n\">{v}</td></tr>" for k, v in items)
    more = f"<caption>{esc(caption)}</caption>" if caption else ""
    return f"<div class=\"tw\"><table>{more}<tbody>{rows}</tbody></table></div>\n"


def paper_page(r, rels, jsonld_node, cite, citation_text):
    k = r["key"]
    title = r.get("title") or r.get("citation_printed") or r["id"]
    au = []
    for a in r.get("authors") or []:
        t = esc(author_text(a))
        if a.get("orcid"):
            t += f' <a href="https://orcid.org/{esc(a["orcid"])}" title="ORCID">(ORCID {esc(a["orcid"])})</a>'
        au.append(t)
    b = [f"<h1>{esc(title)}</h1>\n"]
    if au:
        b.append(f"<p>{', '.join(au)}</p>\n")
    src = " ".join(x for x in [esc(r.get("container_title")), esc(r.get("volume")),
                               f"({esc(r['issue'])})" if r.get("issue") else "",
                               esc(r.get("pages") or r.get("article_number")),
                               f"({r['year']})" if r.get("year") else ""] if x)
    if src:
        b.append(f"<p class=\"lead\">{src}{' · ' + esc(r['publisher']) if r.get('publisher') else ''}</p>\n")
    if r.get("doi"):
        b.append(f"<p>DOI: <a href=\"https://doi.org/{esc(r['doi'])}\">{esc(r['doi'])}</a></p>\n")
    else:
        b.append("<p class=\"muted\">No DOI: id derived from the normalised first author, year and title / printed "
                 f"reference (<code>{esc(r['id'])}</code>).</p>\n")
    if r.get("citation_printed"):
        b.append(f"<p class=\"muted\">Reference as printed in ORCHID's sources: {esc(r['citation_printed'])}</p>\n")
    b.append(coins(r) + "\n")
    b.append(f"<p>Data: <a href=\"../../papers/{k}.json\">JSON record</a> · id <code>{esc(r['id'])}</code> · "
             f"metadata: {esc(r.get('metadata_source'))}</p>\n")
    if r.get("doi_check") and r["doi_check"].get("status") == "mismatch_suspect":
        b.append("<p><strong>Check:</strong> the Crossref record of this DOI differs from the reference printed in "
                 f"ORCHID's source ({esc(', '.join(r['doi_check']['issues']))}); the DOI may be mis-resolved.</p>\n")
    b.append("<h2>Roles in ORCHID</h2>\n<ul>\n")
    for ro in r["roles"]:
        d = r["role_details"].get(ro) or {}
        extra = ""
        if ro == "tdb_source_assessment":
            extra = ": " + ", ".join(esc(t["file"]) + (f" ({esc(sys_label(t['system']))})" if t.get("system") else "")
                                     for t in d.get("tdb_files", []))
        elif ro == "handbook_diagram_source":
            extra = ": " + ", ".join(f"<a href=\"../systems/{esc(s)}.html\">{esc(sys_label(s))}</a>" for s in d.get("systems", []))
        elif ro == "discovery_lead":
            extra = f": rank {esc(d.get('rank'))}, score {esc(d.get('score'))}; systems named in title: " +                     (", ".join(esc(sys_label(s)) for s in d.get("systems_named_in_title") or []) or "none")
        elif ro == "corpus_assessment_read":
            extra = f": {d.get('values_released', 0)} released values"
        b.append(f"<li>{esc(ROLE_LABEL.get(ro, ro))}{extra}</li>\n")
    b.append("</ul>\n")
    if r["systems"]:
        b.append("<h2>Systems</h2>\n<p>" + ", ".join(
            f"<a href=\"../systems/{esc(s)}.html\">{esc(sys_label(s))}</a>" for s in r["systems"]) + "</p>\n")
    es = r["evidence_summary"]
    if es:
        b.append("<h2>Evidence in ORCHID</h2>\n")
        rep = es.get("reported")
        if rep:
            b.append(f"<p>{rep['n_values']} values read from this paper ({rep['n_values_accepted']} accepted), from "
                     f"{rep['n_tables']} tables and {rep['n_text_passages']} text passages.</p>\n")
            b.append(_counter_table("By evidence class", rep["by_evidence_class"]))
            b.append(_counter_table("By property group", rep["by_property_group"]))
        cit = es.get("cited_as_source")
        if cit:
            b.append(f"<p>{cit['n_values']} values in {cit['n_citing_papers']} other papers cite this work as their "
                     "source.</p>\n")
            b.append(_counter_table("Cited values by evidence class", cit["by_evidence_class"]))
        if es.get("espei_datasets"):
            b.append(f"<p>ESPEI datasets derived: {len(es['espei_datasets'])} "
                     f"(<code>{esc(es['espei_datasets'][0])}</code>{' …' if len(es['espei_datasets']) > 1 else ''}).</p>\n")
        if es.get("phase_registry_phases"):
            b.append("<p>Phase-registry phases citing it: " + ", ".join(
                esc(x) for x in es["phase_registry_phases"][:20]) + (" …" if len(es["phase_registry_phases"]) > 20 else "")
                + "</p>\n")
        if es.get("assessed_assemblage_records"):
            b.append(f"<p>Assessed phase-assemblage records citing it: {es['assessed_assemblage_records']}.</p>\n")
    av = r["availability"]
    oa = av.get("open_access_license")
    b.append("<h2>Availability</h2>\n<p>Read by ORCHID: " + esc(av["read_by_orchid"]) +
             (" (" + esc(", ".join(av["read_detail"])) + ")" if av["read_detail"] else "") + ". Open-access licence "
             "(Crossref <code>license</code> field): " +
             (f"<a href=\"{esc(oa['url'])}\">{esc(oa['url'])}</a> ({esc(oa.get('content_version'))})" if oa else "none deposited")
             + ".")
    if r.get("is_referenced_by_count") is not None:
        b.append(f" Crossref is-referenced-by count: {r['is_referenced_by_count']} "
                 f"(retrieved {esc(r['is_referenced_by_count_retrieved'])}).")
    b.append("</p>\n")
    if rels:
        b.append("<h2>Relations</h2>\n")
        groups = {}
        for x in rels:
            groups.setdefault((x["type"], x["direction"]), []).append(x)
        for g in sorted(groups):
            xs = groups[g]
            b.append(f"<h3>{esc(EDGE_LABEL.get(g, g[0] + ' (' + g[1] + ')'))} ({len(xs)})</h3>\n<ul>\n")
            for x in xs[:REL_MAX]:
                o = x["other"]
                if o.get("key"):
                    lab = f"<a href=\"{esc(o['key'])}.html\">{esc(o.get('title') or o['id'])}</a>" +                           (f" ({o['year']})" if o.get("year") else "")
                else:
                    lab = f"<code>{esc(o['id'])}</code>"
                n = f" — {x['count']} values" if x.get("count") else ""
                b.append(f"<li>{lab}{n} <small>[{esc(x['provenance']['source'])}]</small></li>\n")
            if len(xs) > REL_MAX:
                b.append(f"<li>… {len(xs) - REL_MAX} more in the <a href=\"../../papers/{k}.json\">JSON record</a></li>\n")
            b.append("</ul>\n")
    b.append("<aside><h2>Cite</h2><p>Cite the paper itself when you use its data: " + esc(citation_text(r)) +
             "</p><p>Cite this bibliography: " + esc(cite["text"]) + "</p></aside>\n")
    desc = (f"{title}" + (f" ({r['year']})" if r.get("year") else "") + ": roles in ORCHID ("
            + ", ".join(r["roles"]) + "), systems, evidence counts and relations.")
    return page(f"{title} — ORCHID bibliography", desc[:300], f"{DATA_URL}bibliography/pages/papers/{k}.html",
                "".join(b), jsonld_node, depth=2)


def system_page(s, ids, by_id, sysdoc_assessments, cite):
    b = [f"<h1>{esc(sys_label(s))}: papers in ORCHID</h1>\n",
         f"<p>Data: <a href=\"../../systems/{esc(s)}.json\">JSON</a> · <a href=\"../index.html\">all systems</a></p>\n"]
    if sysdoc_assessments:
        b.append("<h2>Assessments and the data papers they used</h2>\n<div class=\"tw\"><table><thead><tr>"
                 "<th scope=\"col\">Assessment</th><th scope=\"col\">Year</th><th scope=\"col\">TDB</th>"
                 "<th scope=\"col\">Data papers used</th><th scope=\"col\">Values from them</th>"
                 "<th scope=\"col\">Reassesses</th></tr></thead><tbody>\n")
        for a in sysdoc_assessments:
            ra = ", ".join(f"<a href=\"../papers/{esc(by_id[x]['key'])}.html\">{esc(by_id[x].get('year') or x)}</a>"
                           for x in a["reassesses"] if x in by_id)
            b.append(f"<tr><td><a href=\"../papers/{esc(a['key'])}.html\">{esc(a.get('title') or a['id'])}</a></td>"
                     f"<td class=\"n\">{esc(a.get('year'))}</td><td>{esc(', '.join(a['tdb_files']))}</td>"
                     f"<td class=\"n\">{a['n_data_papers_used']}</td><td class=\"n\">{a['n_values_from_cited_papers']}</td>"
                     f"<td>{ra}</td></tr>\n")
        b.append("</tbody></table></div>\n")
    for ro, lst in sorted(ids.items()):
        b.append(f"<h2>{esc(ROLE_LABEL.get(ro, ro))} ({len(lst)})</h2>\n<ul>\n")
        for i in sorted(lst, key=lambda i: (by_id[i].get("year") or 0, i)):
            r = by_id[i]
            b.append(f"<li><a href=\"../papers/{esc(r['key'])}.html\">{esc(r.get('title') or r.get('citation_printed') or i)}</a>"
                     + (f" ({r['year']})" if r.get("year") else "") + "</li>\n")
        b.append("</ul>\n")
    n = len({i for v in ids.values() for i in v})
    ld = {"@context": "https://schema.org/", "@type": "Dataset", "name": f"ORCHID bibliography: {sys_label(s)}",
          "description": f"{n} papers about the {sys_label(s)} system known to ORCHID, by role",
          "url": f"{DATA_URL}bibliography/pages/systems/{s}.html",
          "distribution": {"@type": "DataDownload", "contentUrl": f"{DATA_URL}bibliography/systems/{s}.json",
                           "encodingFormat": "application/json"},
          "license": "https://creativecommons.org/licenses/by/4.0/"}
    return page(f"{sys_label(s)} papers — ORCHID bibliography", f"{n} papers about the {sys_label(s)} system known to "
                "ORCHID: assessments, data papers, handbook sources and discovery leads.",
                f"{DATA_URL}bibliography/pages/systems/{s}.html", "".join(b), ld, depth=2)


def index_page(systems, records, cite):
    by_order = {}
    for s in systems:
        by_order.setdefault(s.count("-") + 1, []).append(s)
    b = ["<h1>ORCHID bibliography</h1>\n",
         f"<p class=\"lead\">{len(records)} papers ORCHID knows about, with their roles, evidence and relations; "
         f"{len(systems)} chemical systems.</p>\n",
         "<p>Data: <a href=\"../papers.jsonl\">papers.jsonl</a> · <a href=\"../papers.csv\">CSV</a> · "
         "<a href=\"../papers.parquet\">Parquet</a> · <a href=\"../edges.jsonl\">edges</a> · "
         "<a href=\"../orchid-bibliography.bib\">BibTeX</a> · <a href=\"../orchid-bibliography.csl.json\">CSL-JSON</a> · "
         "<a href=\"../papers.jsonld\">JSON-LD</a> · <a href=\"../index.json\">index.json</a> · "
         "<a href=\"../README.md\">README</a></p>\n"]
    for o in sorted(by_order):
        b.append(f"<h2>{o}-element systems ({len(by_order[o])})</h2>\n<p>" + " · ".join(
            f"<a href=\"systems/{esc(s)}.html\">{esc(sys_label(s))}</a>" for s in by_order[o]) + "</p>\n")
    ld = {"@context": "https://schema.org/", "@type": "Dataset", "name": "ORCHID bibliography",
          "url": f"{DATA_URL}bibliography/pages/index.html", "license": "https://creativecommons.org/licenses/by/4.0/",
          "identifier": f"https://doi.org/{cite['doi']}"}
    return page("ORCHID bibliography", f"{len(records)} papers known to ORCHID with roles, evidence and relations.",
                f"{DATA_URL}bibliography/pages/index.html", "".join(b), ld, depth=1)


def render(records, rel, by_id, systems, cite, jsonld_paper, cites_out, citation_text):
    """{relative path under bibliography/: html} for every paper and system and the index page."""
    out = {}
    for r in records:
        rels = []
        for direction, e in rel.get(r["id"], []):
            other = e["dst"] if direction == "out" else e["src"]
            o = by_id.get(other)
            ob = {"id": other, "key": o["key"], "title": o.get("title") or o.get("citation_printed"),
                  "year": o.get("year")} if o else {"id": other}
            rels.append({"direction": direction, "type": e["type"], "other": ob, "count": e["count"],
                         "provenance": e["provenance"]})
        rels.sort(key=lambda x: (x["type"], x["direction"], -(x["count"] or 0), x["other"]["id"]))
        node = jsonld_paper(r, sorted(cites_out.get(r["id"], [])))
        node["@context"] = "https://schema.org/"
        out[f"pages/papers/{r['key']}.html"] = paper_page(r, rels, node, cite, citation_text)
    for s, ids in systems.items():
        doc = json.loads(_SYSDOCS[s]) if s in _SYSDOCS else None
        out[f"pages/systems/{s}.html"] = system_page(s, ids, by_id, (doc or {}).get("assessments") or [], cite)
    out["pages/index.html"] = index_page(list(systems), records, cite)
    out["pages/style.css"] = STYLE.strip() + "\n"
    return out


_SYSDOCS: dict = {}


def pages_from_files(data_dir):
    """Hook for tools/build_static_pages.py: rebuild the pages from the published bibliography files.
    Returns {"bibliography/pages/...": html}."""
    import collections
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import build_bibliography as bb
    root = os.path.join(data_dir, "bibliography")
    with open(os.path.join(root, "papers.jsonl"), encoding="utf-8") as f:
        records = [json.loads(x) for x in f]
    with open(os.path.join(root, "edges.jsonl"), encoding="utf-8") as f:
        edges = [json.loads(x) for x in f]
    with open(os.path.join(root, "systems_index.json"), encoding="utf-8") as f:
        systems = json.load(f)["systems"]
    with open(os.path.join(root, "index.json"), encoding="utf-8") as f:
        cite = json.load(f)["cite_as"]
    by_id = {r["id"]: r for r in records}
    rel = collections.defaultdict(list)
    cites_out = collections.defaultdict(list)
    for e in edges:
        rel[e["src"]].append(("out", e))
        rel[e["dst"]].append(("in", e))
        if e["type"] in ("cites", "cited_by_crossref") and e["dst"] in by_id:
            u = by_id[e["dst"]]["urls"]["json"]
            if u not in cites_out[e["src"]]:
                cites_out[e["src"]].append(u)
    _SYSDOCS.clear()
    for s in systems:
        with open(os.path.join(root, "systems", f"{s}.json"), encoding="utf-8") as f:
            _SYSDOCS[s] = f.read()
    out = render(records, rel, by_id, systems, cite, bb.jsonld_paper, cites_out, bb.citation_text)
    return {f"bibliography/{k}": v for k, v in out.items()}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="write bibliography/pages/ from the published bibliography files")
    ap.add_argument("--data", default=os.environ.get("ORCHID_DATA") or os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "orchid-data"))
    a = ap.parse_args()
    for rel_path, text in sorted(pages_from_files(a.data).items()):
        p = os.path.join(a.data, rel_path)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            f.write(text)
