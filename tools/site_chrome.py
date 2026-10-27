#!/usr/bin/env python3
"""Shared page chrome of orchid-project.org: one header (navigation + release status), one footer, one stylesheet.

Every website page carries three generated blocks between markers,

  <!-- orchid-chrome-css:begin -->  ... <!-- orchid-chrome-css:end -->   (end of <head>)
  <!-- orchid-header:begin -->      ... <!-- orchid-header:end -->       (start of <body>)
  <!-- orchid-footer:begin -->      ... <!-- orchid-footer:end -->       (end of <body>)

and tools/build_site_numbers.py rewrites them from the definitions below (`--check` fails when a page is stale or a
marker is missing). Change the navigation, the status statement or the footer here, never in a page. The chrome is
self-styled (classes prefixed `o`, colours of its own) so it looks the same on pages with other stylesheets
(the catalog tools). Standard library only.
"""
from __future__ import annotations

import html
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402

# (key, label, target): target is a path from the site root ("" = home) or an absolute URL
NAV = (
    ("explorer", "Data explorer", "data-explorer/"),
    ("catalog", "Catalog", "catalog/"),
    ("diagrams", "Phase diagrams", "phase-diagrams/"),
    ("phases", "Phases", "phases/"),
    ("bibliography", "Bibliography", "https://data.orchid-project.org/bibliography/pages/index.html"),
    ("hulls", "Hulls", "catalog/binary-thermo.html"),
    ("data", "Data", "dataset/"),
    ("downloads", "Downloads", "#data"),
    ("releases", "Releases", "#releases"),
    ("cite", "Cite", "#cite"),
    ("workflow", "Workflow", "workflow/"),
    ("about", "About", "about/"),
    ("contribute", "Contribute", "contribute/"),
)
# page (path from the site root) -> NAV key marked aria-current="page"
CURRENT = {"data-explorer/index.html": "explorer", "catalog/index.html": "catalog", "phase-diagrams/index.html": "diagrams",
           "phases/index.html": "phases", "catalog/binary-thermo.html": "hulls", "dataset/index.html": "data",
           "workflow/index.html": "workflow", "about/index.html": "about", "contribute/index.html": "contribute"}
# every page that carries the chrome (build_site_numbers.py checks them all)
PAGES = ("index.html", "data-explorer/index.html", "catalog/index.html", "phase-diagrams/index.html", "phases/index.html",
         "catalog/binary-thermo.html", "catalog/mixing-energy.html", "catalog/model-audit/index.html",
         "catalog/compatibility/unary-compatibility-network.html", "dataset/index.html", "dft-parity/index.html",
         "workflow/index.html", "about/index.html", "contribute/index.html", "best-practices.html", "404.html")

# Release status, shown under the navigation of every page: orchid_cite.STATUS (two sentences), "release notes" linked
# to the Releases section of the home page
STATUS_LEAD = "Version 1, machine-extracted:"
assert orchid_cite.STATUS.startswith(STATUS_LEAD) and orchid_cite.STATUS.endswith("the release notes.")
STATUS_TEXT = orchid_cite.STATUS[len(STATUS_LEAD):].strip()[:-len("release notes.")] + "{release notes}."
FUNDING = "This work was supported by the U.S. National Science Foundation under awards OAC-2608825 and 2534344."
GITHUB = "https://github.com/orchid-calphad"

CSS = """.ohead,.ofoot{--o-bg:#fbfaf7;--o-fg:#1d1b20;--o-muted:#5e5a66;--o-line:#e2dfe6;--o-card:#fff;--o-accent:#7b3fa0;--o-accent2:#f2eaf7}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]) .ohead,:root:not([data-theme="light"]) .ofoot{--o-bg:#141218;--o-fg:#ece8f1;--o-muted:#a69fb0;--o-line:#2e2a35;--o-card:#1c1922;--o-accent:#c79ae6;--o-accent2:#2a2133}}
:root[data-theme="dark"] .ohead,:root[data-theme="dark"] .ofoot{--o-bg:#141218;--o-fg:#ece8f1;--o-muted:#a69fb0;--o-line:#2e2a35;--o-card:#1c1922;--o-accent:#c79ae6;--o-accent2:#2a2133}
header.ohead{display:block;position:static;margin:0;padding:0;background:var(--o-bg);color:var(--o-fg);border:0;border-bottom:1px solid var(--o-line);box-shadow:none;font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;text-align:left}
footer.ofoot{display:block;margin:48px 0 0;padding:20px 0 36px;background:var(--o-bg);color:var(--o-muted);border:0;border-top:1px solid var(--o-line);font:14px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif;text-align:left;max-width:none}
.ohead *,.ofoot *{box-sizing:border-box}
.oin{max-width:1180px;margin:0 auto;padding:0 16px}
.obar{display:flex;align-items:center;gap:8px 16px;padding:10px 0;flex-wrap:wrap}
.ohead a.obrand{display:inline-flex;align-items:center;gap:8px;margin:0 auto 0 0;font-weight:700;font-size:1.1rem;color:var(--o-fg);text-decoration:none}
.obrand img{display:block;width:26px;height:26px}
nav.onav{display:block;margin:0;padding:0;font:inherit;background:none;border:0}
.onav ul{display:flex;flex-wrap:wrap;gap:2px 14px;list-style:none;margin:0;padding:0}
.onav li{margin:0;padding:0}
.ohead .onav a{display:inline-block;margin:0;padding:4px 0;color:var(--o-accent);text-decoration:none;font-size:.9rem;font-weight:400;white-space:nowrap;border-bottom:2px solid transparent}
.ohead .onav a:hover{border-bottom-color:var(--o-line)}
.ohead .onav a[aria-current="page"]{color:var(--o-fg);font-weight:600;border-bottom-color:var(--o-accent)}
.ohead button.omenu{display:none;align-items:center;gap:6px;margin:0;padding:5px 12px;font:600 .9rem/1.3 system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--o-accent);background:var(--o-card);border:1px solid var(--o-line);border-radius:8px;cursor:pointer}
.ostatus{margin:0;padding:0 0 10px;color:var(--o-muted);font-size:.82rem;line-height:1.45}
.ostatus b{color:var(--o-fg)}
.ostatus a{text-decoration:underline}
.ohead a,.ofoot a{color:var(--o-accent)}
.ofoot p{margin:0 0 8px;max-width:none}
a.oskip{position:absolute;left:-9999px;top:0;z-index:1000;padding:8px 12px;background:var(--o-card);color:var(--o-accent)}
a.oskip:focus{left:8px;top:8px}
:focus-visible{outline:2px solid #7b3fa0;outline-offset:2px}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]) :focus-visible{outline-color:#c79ae6}}
@media (max-width:1060px){.ohead.ojs button.omenu{display:inline-flex}.ohead.ojs nav.onav{display:none;flex-basis:100%}
.ohead.ojs.oopen nav.onav{display:block}.ohead.ojs.oopen .onav ul{flex-direction:column;gap:0;padding:4px 0 6px}
.ohead.ojs .onav a{padding:7px 0;font-size:1rem}}
@media (max-width:600px){.ostatus{font-size:.78rem}}"""

SCRIPT = ("(function(){var s=document.currentScript,h=s&&s.previousElementSibling;if(!h||!h.classList.contains('ohead'))return;"
          "var b=h.querySelector('.omenu');h.classList.add('ojs');b.addEventListener('click',function(){"
          "var o=h.classList.toggle('oopen');b.setAttribute('aria-expanded',o?'true':'false');});})();")

BLOCK = re.compile(r"<!-- orchid-(chrome-css|header|footer):begin[^>]*-->.*?<!-- orchid-\1:end -->", re.S)
NOTE = " (tools/site_chrome.py; rewritten by tools/build_site_numbers.py)"


def prefix(page: str) -> str:
    """Relative path from `page` to the site root ('/' for the 404 page, which is served at any depth)."""
    if page == "404.html":
        return "/"
    depth = page.count("/")
    return "../" * depth if depth else "./"


def href(target: str, page: str) -> str:
    if target.startswith(("http://", "https://")):
        return target
    if target.startswith("#"):
        return target if page == "index.html" else prefix(page) + target
    return prefix(page) + target


def css_block() -> str:
    return "<!-- orchid-chrome-css:begin%s -->\n<style id=\"orchid-chrome\">\n%s\n</style>\n<!-- orchid-chrome-css:end -->" % (
        NOTE, CSS)


def status_html(page: str) -> str:
    link = '<a href="%s">release notes</a>' % html.escape(href("#releases", page))
    return '<p class="ostatus"><b>%s</b> %s</p>' % (html.escape(STATUS_LEAD), html.escape(STATUS_TEXT).replace(
        "{release notes}", link))


def header_block(page: str) -> str:
    cur = CURRENT.get(page)
    lis = []
    for key, label, target in NAV:
        attr = ' aria-current="page"' if key == cur else ""
        lis.append('<li><a href="%s"%s>%s</a></li>' % (html.escape(href(target, page)), attr, html.escape(label)))
    home = prefix(page)
    return "\n".join([
        "<!-- orchid-header:begin%s -->" % NOTE,
        '<a class="oskip" href="#ocontent">Skip to content</a>',
        '<header class="ohead"><div class="oin">',
        '<div class="obar"><a class="obrand" href="%s"%s><img src="%sassets/logo/orchid-icon.svg" alt="" width="26" '
        'height="26"><span>ORCHID</span></a>' % (html.escape(home), ' aria-current="page"' if page == "index.html" else "",
                                                  html.escape(home)),
        '<button class="omenu" type="button" aria-expanded="false" aria-controls="onav">Menu</button>',
        '<nav class="onav" id="onav" aria-label="Site"><ul>',
        "\n".join(lis),
        "</ul></nav></div>",
        status_html(page),
        "</div></header>",
        "<script>%s</script>" % SCRIPT,
        '<span id="ocontent"></span>',
        "<!-- orchid-header:end -->"])


def footer_block(page: str) -> str:
    p = prefix(page)
    links = [("About", p + "about/"), ("Credits: software and data citations", p + "about/#software-data-credits"),
             ("CITATIONS.md", p + "CITATIONS.md"), ("Contribute", p + "contribute/"), ("llms.txt", p + "llms.txt"),
             ("Sitemap", p + "sitemap.xml"), ("GitHub", GITHUB)]
    return "\n".join([
        "<!-- orchid-footer:begin%s -->" % NOTE,
        '<footer class="ofoot"><div class="oin">',
        "<p>%s</p>" % html.escape(FUNDING),
        "<p>ORCHID is developed at Texas A&amp;M University with Citrine Informatics and Oak Ridge National Laboratory. "
        "Site content: MIT; ORCHID data: CC BY 4.0; database files remain the property of their authors "
        '(<a href="%sdata/NOTICE.md">NOTICE</a>).</p>' % html.escape(p),
        "<p>%s</p>" % " · ".join('<a href="%s">%s</a>' % (html.escape(u), html.escape(t)) for t, u in links),
        "</div></footer>",
        "<!-- orchid-footer:end -->"])


def blocks(page: str) -> dict:
    return {"chrome-css": css_block(), "header": header_block(page), "footer": footer_block(page)}


def apply(page: str, text: str):
    """(new text, missing markers): every chrome block rewritten from the definitions above."""
    want = blocks(page)
    seen = set()

    def sub(m):
        seen.add(m.group(1))
        return want[m.group(1)]
    new = BLOCK.sub(sub, text)
    return new, sorted(set(want) - seen)
