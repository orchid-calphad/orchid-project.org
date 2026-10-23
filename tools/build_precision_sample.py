#!/usr/bin/env python3
"""Precision-sample review kit for the ORCHID literature data: draw once, then freeze and verify.

The sample is FIXED. It was drawn once (seed 20261022, 300 items) from one literature export and distributed to
reviewers (index.html and the packets); it is never regenerated from a newer export. Files in
review/precision-sample/:

  sample.json    the sample (records + strata + allocation), canonical JSON; frozen
  index.html     the self-contained review page embedding sample.json; frozen
  packets/       the reviewer packets cut from the page; frozen
  frozen.json    the freeze record: SHA-256 of the files above, the population snapshot the sample was drawn from
                 (export SHA-256, papers SHA-256, release, accepted population, orchid-data commit), and every sampled
                 record id with the field values the reviewers see (per review question) and their hash

Commands (standard library only):

  build_precision_sample.py --check [--drift-json PATH]
      verify (a) sample.json, index.html and the packets are byte-identical to the freeze record, and (b) every
      sampled record id still exists in the current export (data/literature/orchid-literature.jsonl.gz). Records whose
      released values changed since sampling are listed (and written to PATH) for the scorer; that is not a failure.
      Running the script with no command does the same.
  build_precision_sample.py freeze [--verify-draw]
      write frozen.json for an existing sample (refuses to change an existing freeze record). --verify-draw re-draws
      the sample from the population snapshot in the orchid-data git history and requires identical bytes.
  build_precision_sample.py draw --out NEW_DIR [--seed S] [-n N]
      draw a NEW sample from the current export into a new directory (refuses an existing sample there) and freeze it.

Strata: evidence class x record kind (invariant_table / general_table / running_text). The sample size
(300) is allocated proportionally to stratum size with a minimum of 15 per non-empty stratum (all of a
stratum when it holds fewer than 15). Records are sorted by id within a stratum and drawn with
random.Random(seed).sample, so the same input file and seed give the same bytes.

Results downloaded from the page are scored by tools/score_precision_sample.py against the frozen sample.
"""
import argparse
import gzip
import hashlib
import json
import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INPUT = ROOT / "data" / "literature" / "orchid-literature.jsonl.gz"
PAPERS = ROOT / "data" / "literature" / "papers.json"
OUTDIR = ROOT / "review" / "precision-sample"
SEED = 20261022
N = 300
MIN_PER_STRATUM = 15
KINDS = ("invariant_table", "general_table", "running_text")
EV_ORDER = ("experimental", "model", "dft", "unspecified")
SCHEMA = "orchid-precision-sample/1"
FREEZE_SCHEMA = "orchid-precision-sample-freeze/1"
DATA_REPO = ROOT.parent / "orchid-data"
# what a reviewer judges, per review question (the page shows these fields); compared with the current release
REVIEWED = {
    "value": ("system", "property", "value", "value_range", "value_expression", "value_printed", "phase_compositions",
              "temperature_K", "temperature_printed"),
    "unit": ("unit", "unit_printed", "unit_basis"),
    "evidence_class": ("evidence_class",),
    "phases": ("phases",),
    "citation": ("paper", "paper_doi", "cited_sources"),
}


def canonical(obj):
    """Sorted keys, 2-space indent, UTF-8, trailing newline (same bytes as the page's canonicalJSON)."""
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"


def allocate(pop, n, minimum):
    """Proportional allocation (largest remainder) with a per-stratum minimum, capped at stratum size."""
    fixed = {}
    while True:
        free = {k: v for k, v in pop.items() if k not in fixed}
        left = n - sum(fixed.values())
        total = sum(free.values())
        quota = {k: left * v / total for k, v in free.items()}
        floor = {k: min(minimum, pop[k]) for k in free}
        low = {k for k in free if quota[k] < floor[k] or quota[k] > pop[k]}
        if not low:
            break
        for k in low:
            fixed[k] = floor[k] if quota[k] < floor[k] else pop[k]
    base = {k: math.floor(q) for k, q in quota.items()}
    rest = left - sum(base.values())
    for k in sorted(quota, key=lambda k: (-(quota[k] - base[k]), k))[:rest]:
        base[k] += 1
    base.update(fixed)
    return base


def stratum_key(r):
    return "%s|%s" % (r["evidence_class"], r["record_kind"])


def build(input_path, papers_path, seed, n):
    raw = input_path.read_bytes()
    pop_records = {}
    for line in gzip.decompress(raw).decode("utf-8").splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        if not (r.get("flags") or {}).get("extraction_accepted"):
            continue
        if r.get("record_kind") not in KINDS:
            continue
        pop_records.setdefault(stratum_key(r), []).append(r)
    pop = {k: len(v) for k, v in pop_records.items()}
    alloc = allocate(pop, n, MIN_PER_STRATUM)
    papers = json.loads(papers_path.read_text(encoding="utf-8"))

    def paper_info(pid):
        p = papers.get(pid) or {}
        return {"id": pid, "citation": p.get("citation"), "doi": p.get("doi"), "doi_url": p.get("doi_url"),
                "title": p.get("title"), "year": p.get("year")}

    rng = random.Random(seed)
    items = []
    order = sorted(pop, key=lambda k: (EV_ORDER.index(k.split("|")[0]) if k.split("|")[0] in EV_ORDER else 99,
                                       KINDS.index(k.split("|")[1])))
    for k in order:
        recs = sorted(pop_records[k], key=lambda r: r["id"])
        for r in rng.sample(recs, alloc[k]):
            items.append({"stratum": k, "record": r, "paper": paper_info(r["paper"])})
    for i, it in enumerate(items, 1):
        it["n"] = i
    strata = [{"stratum": k, "evidence_class": k.split("|")[0], "record_kind": k.split("|")[1],
               "population": pop[k], "sampled": alloc[k]} for k in order]
    release = None
    ix = input_path.parent / "index.json"
    if ix.exists():
        try:
            release = json.loads(ix.read_text(encoding="utf-8")).get("release")
        except ValueError:
            release = None
    return {
        "schema": SCHEMA,
        "description": "Stratified random sample of accepted ORCHID literature values for a human precision check.",
        "input": {"file": "data/literature/orchid-literature.jsonl.gz", "sha256": hashlib.sha256(raw).hexdigest(),
                  "papers_sha256": hashlib.sha256(papers_path.read_bytes()).hexdigest(), "release": release,
                  "population_accepted": sum(pop.values())},
        "design": {"seed": seed, "n": n, "min_per_stratum": MIN_PER_STRATUM,
                   "strata": "evidence_class x record_kind", "allocation": "proportional, largest remainder, minimum per non-empty stratum",
                   "draw": "records sorted by id within stratum; random.Random(seed).sample, strata in fixed order"},
        "strata": strata,
        "items": items,
    }


PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="robots" content="noindex, nofollow">
<title>Precision sample review</title>
<link rel="icon" type="image/svg+xml" href="../../assets/logo/orchid-icon.svg">
<style>
:root{--bg:#fbfaf7;--fg:#1d1b20;--muted:#5e5a66;--line:#e2dfe6;--card:#fff;--accent:#7b3fa0;--accent2:#f2eaf7;
--exp:#1d7448;--exp-bg:#e3f3ea;--mod:#2c58a8;--mod-bg:#e5ecf9;--dft:#a3520a;--dft-bg:#fbeedf;--uns:#5f5a66;--uns-bg:#ecebee;--warn:#9a3412}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#141218;--fg:#ece8f1;--muted:#a69fb0;--line:#2e2a35;--card:#1c1922;--accent:#c79ae6;--accent2:#2a2133;
--exp:#74d3a2;--exp-bg:#183326;--mod:#9dbcf3;--mod-bg:#1b2840;--dft:#f2ad6b;--dft-bg:#3a2612;--uns:#bdb7c6;--uns-bg:#2a2730;--warn:#fdba74}}
:root[data-theme="dark"]{--bg:#141218;--fg:#ece8f1;--muted:#a69fb0;--line:#2e2a35;--card:#1c1922;--accent:#c79ae6;--accent2:#2a2133;
--exp:#74d3a2;--exp-bg:#183326;--mod:#9dbcf3;--mod-bg:#1b2840;--dft:#f2ad6b;--dft-bg:#3a2612;--uns:#bdb7c6;--uns-bg:#2a2730;--warn:#fdba74}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:var(--accent)}
.wrap{max-width:1040px;margin:0 auto;padding:0 16px}
header{border-bottom:1px solid var(--line)}
nav{display:flex;flex-wrap:wrap;gap:4px 18px;align-items:center;padding:14px 0}
nav .brand{font-weight:700;font-size:1.15rem;margin-right:auto;text-decoration:none;color:var(--fg);display:inline-flex;align-items:center;gap:8px}
h1{font-size:clamp(1.4rem,4vw,1.9rem);margin:24px 0 6px;line-height:1.2}
.lead{color:var(--muted);max-width:800px;margin:0 0 6px}
.note{background:var(--accent2);border-radius:8px;padding:10px 12px;margin:12px 0;font-size:.95rem}
.bar{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:14px 0}
.bar label{font-size:.9rem}
input[type=text],select,textarea{font:inherit;color:var(--fg);background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px 8px;max-width:100%}
textarea{width:100%;min-height:64px}
button{font:inherit;color:var(--fg);background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px 12px;cursor:pointer}
button.primary{background:var(--accent);color:var(--bg);border-color:var(--accent)}
button:disabled{opacity:.45;cursor:default}
.progress{font-variant-numeric:tabular-nums;font-weight:600}
.meter{height:6px;background:var(--line);border-radius:3px;overflow:hidden;flex:1 1 160px;min-width:120px}
.meter i{display:block;height:100%;background:var(--accent);width:0}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:12px 0}
dl{display:grid;grid-template-columns:minmax(110px,170px) 1fr;gap:4px 14px;margin:0}
dt{color:var(--muted);font-size:.88rem}
dd{margin:0;overflow-wrap:anywhere}
.value{font-size:1.25rem;font-weight:600;font-variant-numeric:tabular-nums}
blockquote{margin:0;padding:8px 12px;border-left:3px solid var(--accent);background:var(--accent2);border-radius:0 6px 6px 0}
.badge{display:inline-block;font-size:.75rem;font-weight:600;padding:1px 8px;border-radius:999px}
.b-experimental{background:var(--exp-bg);color:var(--exp)}.b-model{background:var(--mod-bg);color:var(--mod)}
.b-dft{background:var(--dft-bg);color:var(--dft)}.b-unspecified{background:var(--uns-bg);color:var(--uns)}
fieldset{border:1px solid var(--line);border-radius:8px;margin:8px 0;padding:6px 10px}
legend{font-size:.9rem;font-weight:600;padding:0 4px}
fieldset label{display:inline-flex;align-items:center;gap:4px;margin:2px 14px 2px 0;font-size:.95rem}
.q{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:4px 12px}
.muted{color:var(--muted)}
.warn{color:var(--warn)}
code{font-size:.85em;background:var(--accent2);padding:1px 5px;border-radius:4px;overflow-wrap:anywhere}
details pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;max-height:340px;overflow:auto;background:var(--accent2);padding:8px;border-radius:6px}
footer{border-top:1px solid var(--line);margin-top:40px;padding:18px 0 36px;color:var(--muted);font-size:.88rem}
</style>
</head>
<body>
<header><div class="wrap"><nav>
  <a class="brand" href="../../"><img src="../../assets/logo/orchid-icon.svg" alt="" width="28" height="28">ORCHID</a>
  <span class="muted">Internal review</span>
</nav></div></header>
<main class="wrap">
<h1>Precision sample: check accepted values against their papers</h1>
<p class="lead">A stratified random sample of <b id="n-items">…</b> accepted literature values. For each one, open the paper, find the
value at the given page, table and row (or the quoted sentence) and answer the questions. Answers are saved in this browser as you go;
use <b>Download results</b> to send them back.</p>
<p class="note"><b>To verify, open the paper through its DOI link</b> (publisher page, institutional access). This page holds no page
images or PDFs. Mark <i>can't tell</i> when the paper is not accessible or the item cannot be judged; <i>n/a</i> when the field is empty
and should be.</p>

<div class="bar">
  <label for="reviewer">Reviewer name</label>
  <input type="text" id="reviewer" autocomplete="name" placeholder="Your name" style="flex:1 1 200px">
</div>
<div class="bar">
  <span class="progress" id="progress">0 / 0 answered</span>
  <span class="meter" aria-hidden="true"><i id="meter"></i></span>
</div>
<div class="bar">
  <button type="button" id="prev">← Previous</button>
  <label for="jump" class="muted">Item</label>
  <select id="jump" aria-label="Jump to item" style="flex:1 1 180px;min-width:0"></select>
  <button type="button" id="next">Next →</button>
  <button type="button" id="next-open">Next unanswered</button>
</div>

<article class="card" id="item" aria-live="polite"></article>

<section class="card" id="answers" aria-label="Your answers">
  <div class="q" id="questions"></div>
  <label for="note" style="display:block;margin-top:8px;font-size:.9rem;font-weight:600">Note (what is wrong, where you looked)</label>
  <textarea id="note"></textarea>
  <p class="muted" id="saved" style="font-size:.85rem;margin:6px 0 0"></p>
</section>

<div class="bar">
  <button type="button" class="primary" id="download">Download results</button>
  <button type="button" id="clear">Clear my answers</button>
  <span class="muted" id="status"></span>
</div>
<p class="muted" style="font-size:.85rem">Sample <code id="sha"></code> · seed <span id="seed"></span> · strata: evidence class × source kind.
Results are scored by <code>tools/score_precision_sample.py</code>.</p>
</main>
<footer><div class="wrap">Internal review page; not linked from the site. Quotes are single sentences reproduced under fair use so that values can be checked; source papers are not redistributed.</div></footer>
<script type="application/json" id="sample-data">__SAMPLE__</script>
<script>
"use strict";
const SAMPLE_SHA256 = "__SHA__";
const $ = id => document.getElementById(id);
function h(t, a, ...k) { const e = document.createElement(t); if (a) for (const x in a) { if (x === "class") e.className = a[x]; else e.setAttribute(x, a[x]); } k.forEach(c => { if (c != null) e.append(c instanceof Node ? c : String(c)); }); return e; }
const S = JSON.parse($("sample-data").textContent);
const ITEMS = S.items;
const KEY = "orchid-precision-sample:" + SAMPLE_SHA256.slice(0, 16);
const QUESTIONS = [
  {id: "value", label: "Value correct?", opts: ["yes", "no", "cant_tell"]},
  {id: "unit", label: "Unit / basis correct?", opts: ["yes", "no", "cant_tell", "na"]},
  {id: "evidence_class", label: "Evidence class correct?", opts: ["yes", "no", "cant_tell"]},
  {id: "phases", label: "Phases correct?", opts: ["yes", "no", "cant_tell", "na"]},
  {id: "citation", label: "Citation correct?", opts: ["yes", "no", "cant_tell"]},
];
const OPTL = {yes: "yes", no: "no", cant_tell: "can't tell", na: "n/a"};
const KINDL = {invariant_table: "invariant-reaction table", general_table: "general table", running_text: "running text"};
let state = {reviewer: "", answers: {}};
let cur = 0;

function load() {
  try { const s = localStorage.getItem(KEY); if (s) { const v = JSON.parse(s); if (v && typeof v === "object" && v.answers) state = {reviewer: String(v.reviewer || ""), answers: v.answers}; } } catch (e) { /* storage unavailable */ }
  try { const c = Number(localStorage.getItem(KEY + ":cur")); if (c >= 0 && c < ITEMS.length) cur = c; } catch (e) {}
}
function save() {
  try { localStorage.setItem(KEY, JSON.stringify(state)); localStorage.setItem(KEY + ":cur", String(cur)); $("saved").textContent = "Saved in this browser."; }
  catch (e) { $("saved").textContent = "Not saved: browser storage is unavailable. Download results before closing."; }
}
function fmtNum(v) { return typeof v === "number" ? String(v) : (v == null ? "" : String(v)); }
function isAnswered(id) { const a = state.answers[id]; return !!a && QUESTIONS.every(q => a[q.id]); }

function valueText(r) {
  if (r.value != null) return fmtNum(r.value) + (r.unit ? " " + r.unit : "");
  if (r.value_range) return fmtNum(r.value_range.min) + " – " + fmtNum(r.value_range.max) + (r.unit ? " " + r.unit : "");
  if (r.value_expression) return String(r.value_expression);
  if (r.phase_compositions && r.phase_compositions.length) return "phase compositions (below)";
  return "(no scalar value; see phases)";
}
function compText(c) { return Object.keys(c || {}).sort().map(k => k + " " + c[k]).join(", "); }
function badge(ev) { return h("span", {class: "badge b-" + ev}, ev === "dft" ? "DFT" : ev); }
function row(dl, k, v) { if (v == null || v === "" || (Array.isArray(v) && !v.length)) return; dl.append(h("dt", null, k), h("dd", null, v)); }
function safeURL(u) { return /^https:\/\/(dx\.)?doi\.org\//.test(u || "") ? u : null; }

function renderItem() {
  const it = ITEMS[cur], r = it.record, p = it.paper || {}, box = $("item");
  box.replaceChildren();
  box.append(h("p", {class: "muted", style: "margin:0 0 6px;font-size:.9rem"}, "Item " + it.n + " of " + ITEMS.length + " · stratum ", h("b", null, r.evidence_class + " × " + (KINDL[r.record_kind] || r.record_kind)), " · ", h("code", null, r.id)));
  const dl = h("dl");
  row(dl, "System", r.system);
  row(dl, "Property", (r.property || "") + (r.property_printed ? " (printed: " + r.property_printed + ")" : "") + (r.property_group ? " · " + r.property_group : ""));
  row(dl, "Value", h("span", {class: "value"}, valueText(r)));
  if (r.value_printed != null || r.unit_printed) row(dl, "As printed", fmtNum(r.value_printed) + (r.unit_printed ? " " + r.unit_printed : ""));
  if (r.uncertainty && r.uncertainty.basis === "printed") row(dl, "Uncertainty", "± " + fmtNum(r.uncertainty.value) + " (printed)");
  const T = r.temperature_K != null ? fmtNum(r.temperature_K) + " K" + (r.temperature_printed ? " (printed: " + r.temperature_printed + ")" : "") : (r.temperature_printed || "");
  row(dl, "Temperature", T + (r.temperature_basis && r.temperature_basis !== "printed" ? " · basis: " + r.temperature_basis : ""));
  if (r.reaction) row(dl, "Reaction", (r.reaction.printed || "") + (r.reaction.type ? " · " + r.reaction.type : ""));
  if (r.phases && r.phases.length) row(dl, "Phases", r.phases.map(x => x.label + (x.model_phase ? " → " + x.model_phase : "")).join("; "));
  if (r.phase_compositions && r.phase_compositions.length) row(dl, "Phase compositions", h("span", null, ...r.phase_compositions.map((pc, i) => h("span", {style: "display:block"}, (pc.phase || "?") + ": " + compText(pc.composition) + (pc.basis ? " (" + pc.basis + ")" : "")))));
  if (r.overall_composition) row(dl, "Overall composition", compText(r.overall_composition));
  if (r.method) row(dl, "Method", r.method);
  row(dl, "Evidence class", h("span", null, badge(r.evidence_class), r.origin && r.origin.basis ? h("span", {class: "muted", style: "font-size:.85rem"}, " basis: " + r.origin.basis) : null));
  const cit = h("span", null, p.citation || p.id || r.paper);
  const u = safeURL(p.doi_url);
  if (u) cit.append(" ", h("a", {href: u, target: "_blank", rel: "noopener noreferrer"}, "doi:" + p.doi));
  else cit.append(" ", h("span", {class: "warn"}, "(no DOI on record: find the paper by its citation)"));
  row(dl, "Paper", cit);
  if (r.cited_sources && r.cited_sources.length) row(dl, "Cited source(s) of the value", h("span", null, ...r.cited_sources.map(c => { const s = h("span", {style: "display:block"}, (c.label ? "[" + c.label + "] " : "") + (c.reference || c.paper || "")); if (c.doi) s.append(" ", h("a", {href: "https://doi.org/" + encodeURI(c.doi), target: "_blank", rel: "noopener noreferrer"}, "doi")); return s; })));
  const L = r.locator || {};
  if (L.artifact === "table") row(dl, "Where", "page " + fmtNum(L.page) + " · Table " + fmtNum(L.label) + (L.row != null ? " · row " + fmtNum(L.row) : "") + (L.caption ? " — " + L.caption : ""));
  else row(dl, "Where", "page " + fmtNum(r.quote_page != null ? r.quote_page : L.page) + (L.label ? " · " + L.label : ""));
  if (r.quote) row(dl, "Quote", h("blockquote", null, r.quote));
  const fl = r.flags || {}, notes = [];
  if (fl.sign_pdf_page_verified) notes.push("sign corrected from the PDF page");
  if (fl.unsigned_in_print) notes.push("printed without a minus sign; kept as printed");
  if (fl.has_dropped_fields) notes.push("fields not printed in the quote were removed");
  if (notes.length) row(dl, "Flags", notes.join("; "));
  box.append(dl);
  box.append(h("details", {style: "margin-top:10px"}, h("summary", {class: "muted", style: "font-size:.88rem"}, "Full record"), h("pre", null, JSON.stringify(r, null, 2))));
  renderAnswers();
  $("jump").value = String(cur);
  $("prev").disabled = cur === 0; $("next").disabled = cur === ITEMS.length - 1;
}

function renderAnswers() {
  const id = ITEMS[cur].record.id, a = state.answers[id] || {}, qs = $("questions");
  qs.replaceChildren();
  QUESTIONS.forEach(q => {
    const fs = h("fieldset", null, h("legend", null, q.label));
    q.opts.forEach(o => {
      const inp = h("input", {type: "radio", name: "q-" + q.id, value: o, id: "q-" + q.id + "-" + o});
      if (a[q.id] === o) inp.checked = true;
      inp.addEventListener("change", () => setAnswer(q.id, o));
      fs.append(h("label", {for: inp.id}, inp, OPTL[o]));
    });
    qs.append(fs);
  });
  $("note").value = a.note || "";
  $("saved").textContent = "";
}
function ensure(id) { return state.answers[id] || (state.answers[id] = {}); }
function setAnswer(field, v) { const a = ensure(ITEMS[cur].record.id); a[field] = v; a.updated = new Date().toISOString(); save(); progress(); }
function progress() {
  const done = ITEMS.filter(it => isAnswered(it.record.id)).length;
  $("progress").textContent = done + " / " + ITEMS.length + " answered";
  $("meter").style.width = (100 * done / ITEMS.length).toFixed(1) + "%";
  [...$("jump").options].forEach((o, i) => { o.textContent = (isAnswered(ITEMS[i].record.id) ? "✓ " : "· ") + ITEMS[i].n + " · " + ITEMS[i].record.system + " · " + ITEMS[i].record.property; });
}
function go(i) { cur = Math.max(0, Math.min(ITEMS.length - 1, i)); save(); renderItem(); }

function canonicalize(v) {
  if (Array.isArray(v)) return v.map(canonicalize);
  if (v && typeof v === "object") { const o = {}; Object.keys(v).sort().forEach(k => { if (v[k] !== undefined) o[k] = canonicalize(v[k]); }); return o; }
  return v;
}
function canonicalJSON(v) { return JSON.stringify(canonicalize(v), null, 2) + "\n"; }

function results() {
  const answers = [];
  ITEMS.forEach(it => {
    const a = state.answers[it.record.id];
    if (!a) return;
    const o = {id: it.record.id, n: it.n, stratum: it.stratum};
    QUESTIONS.forEach(q => { o[q.id] = a[q.id] || null; });
    o.note = (a.note || "").trim();
    o.updated = a.updated || null;
    answers.push(o);
  });
  return {schema: "orchid-precision-review/1", sample_sha256: SAMPLE_SHA256, sample_schema: S.schema, seed: S.design.seed,
          reviewer: state.reviewer.trim(), timestamp: new Date().toISOString(), n_items: ITEMS.length,
          n_answered_complete: ITEMS.filter(it => isAnswered(it.record.id)).length, answers};
}
function download() {
  if (!state.reviewer.trim()) { $("status").textContent = "Enter your name first."; $("reviewer").focus(); return; }
  const blob = new Blob([canonicalJSON(results())], {type: "application/json"});
  const a = h("a", {href: URL.createObjectURL(blob), download: "precision-review-" + state.reviewer.trim().replace(/[^A-Za-z0-9_-]+/g, "_") + "-" + SAMPLE_SHA256.slice(0, 8) + ".json"});
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  $("status").textContent = "Downloaded.";
}

load();
$("n-items").textContent = ITEMS.length;
$("sha").textContent = SAMPLE_SHA256.slice(0, 16);
$("seed").textContent = S.design.seed;
ITEMS.forEach((it, i) => $("jump").append(h("option", {value: String(i)}, String(it.n))));
$("reviewer").value = state.reviewer;
$("reviewer").addEventListener("input", e => { state.reviewer = e.target.value; save(); });
$("note").addEventListener("input", e => { const a = ensure(ITEMS[cur].record.id); a.note = e.target.value; a.updated = new Date().toISOString(); save(); });
$("prev").addEventListener("click", () => go(cur - 1));
$("next").addEventListener("click", () => go(cur + 1));
$("next-open").addEventListener("click", () => { for (let k = 1; k <= ITEMS.length; k++) { const i = (cur + k) % ITEMS.length; if (!isAnswered(ITEMS[i].record.id)) return go(i); } $("status").textContent = "All items answered."; });
$("jump").addEventListener("change", e => go(Number(e.target.value)));
$("download").addEventListener("click", download);
$("clear").addEventListener("click", () => { if (!confirm("Clear all answers saved in this browser?")) return; state.answers = {}; save(); renderItem(); progress(); });
renderItem(); progress();
</script>
</body>
</html>
"""


def render(sample_text):
    sha = hashlib.sha256(sample_text.encode("utf-8")).hexdigest()
    embedded = sample_text.replace("</", "<\\/").replace("<!--", "<\\!--")
    return PAGE.replace("__SAMPLE__", embedded).replace("__SHA__", sha), sha


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(p):
    return sha256_bytes(Path(p).read_bytes())


def reviewed_view(record, paper_doi):
    """The reviewed fields of a record, grouped by review question (absent fields are null)."""
    r = dict(record, paper_doi=paper_doi)
    out = {}
    for q, fields in REVIEWED.items():
        v = {}
        for f in fields:
            x = r.get(f)
            if f == "cited_sources" and x:
                x = [{k: c.get(k) for k in ("label", "reference", "doi", "paper")} for c in x]
            v[f] = x
        out[q] = v
    return out


def view_sha(view):
    return sha256_bytes(json.dumps(view, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))


def freeze_record(outdir, sample_bytes, population_commit=None, draw_verified=None):
    sample = json.loads(sample_bytes)
    page = outdir / "index.html"
    pk = outdir / "packets"
    files = {"sample.json": sha256_bytes(sample_bytes), "index.html": sha256_file(page)}
    if pk.is_dir():
        for p in sorted(pk.iterdir()):
            if p.is_file():
                files["packets/" + p.name] = sha256_file(p)
    items = []
    for it in sample["items"]:
        view = reviewed_view(it["record"], (it.get("paper") or {}).get("doi"))
        items.append({"n": it["n"], "id": it["record"]["id"], "stratum": it["stratum"], "reviewed": view,
                      "reviewed_sha256": view_sha(view)})
    inp = sample["input"]
    return {
        "schema": FREEZE_SCHEMA,
        "note": "The precision sample is frozen: drawn once and distributed to reviewers; never regenerated. "
                "tools/build_precision_sample.py --check verifies these hashes and that every sampled id is still released.",
        "sample_sha256": files["sample.json"],
        "files_sha256": files,
        "design": sample["design"],
        "population": {"file": inp["file"], "sha256": inp["sha256"], "papers_sha256": inp["papers_sha256"],
                       "release": inp["release"], "population_accepted": inp["population_accepted"],
                       "orchid_data_commit": population_commit, "draw_reproduced": draw_verified},
        "reviewed_fields": {q: list(f) for q, f in REVIEWED.items()},
        "items": items,
    }


def load_release(input_path, papers_path):
    """({record id: record}, papers) of the current literature export; ids merged into another record by the
    de-duplication (`merged_ids`) map to the surviving record under the key ("merged", id)."""
    recs = {}
    for line in gzip.decompress(Path(input_path).read_bytes()).decode("utf-8").splitlines():
        if line.strip():
            r = json.loads(line)
            recs[r["id"]] = r
            for m in r.get("merged_ids") or []:
                recs.setdefault(("merged", m), r)
    papers = json.loads(Path(papers_path).read_text(encoding="utf-8"))
    return recs, papers


def drift(frozen, recs, papers):
    """Per sampled item: status 'unchanged' | 'changed' (with the review questions / fields whose released value
    differs from the reviewed one) | 'missing' (id no longer released). `merged_into` names the surviving record when
    the de-duplication merged the sampled id into another (its values are compared)."""
    out = []
    for it in frozen["items"]:
        r = recs.get(it["id"])
        base = {"n": it["n"], "id": it["id"]}
        if r is None:
            r = recs.get(("merged", it["id"]))
            if r is None:
                out.append(dict(base, status="missing"))
                continue
            base["merged_into"] = r["id"]
        cur = reviewed_view(r, (papers.get(r.get("paper")) or {}).get("doi"))
        if view_sha(cur) == it["reviewed_sha256"]:
            out.append(dict(base, status="unchanged"))
            continue
        ch = {}
        for q, fields in it["reviewed"].items():
            fs = [f for f in fields if fields[f] != cur[q].get(f)]
            if fs:
                ch[q] = {f: {"reviewed": fields[f], "released": cur[q].get(f)} for f in fs}
        out.append(dict(base, status="changed" if ch else "unchanged", questions=sorted(ch), changes=ch))
    return out


def verify(outdir, input_path, papers_path):
    """-> (problems [str], drift list). Problems: frozen files changed, sampled ids missing from the release."""
    fp = outdir / "frozen.json"
    if not fp.exists():
        return ["%s missing (run `freeze` once for the distributed sample)" % fp], []
    frozen = json.loads(fp.read_text(encoding="utf-8"))
    bad = []
    for rel, sha in frozen["files_sha256"].items():
        p = outdir / rel
        if not p.exists():
            bad.append("frozen file missing: " + rel)
        elif sha256_file(p) != sha:
            bad.append("frozen file changed: %s (sha256 %s, frozen %s)" % (rel, sha256_file(p)[:16], sha[:16]))
    pk = outdir / "packets"
    if pk.is_dir():
        for p in sorted(pk.iterdir()):
            if p.is_file() and "packets/" + p.name not in frozen["files_sha256"]:
                bad.append("file not in the freeze record: packets/" + p.name)
    recs, papers = load_release(input_path, papers_path)
    d = drift(frozen, recs, papers)
    bad += ["sampled record no longer released: n=%d %s" % (x["n"], x["id"]) for x in d if x["status"] == "missing"]
    return bad, d


def population_commit(sample, repo=DATA_REPO):
    """orchid-data commit whose literature export has the sample's population SHA-256 (git history, local)."""
    rel = "literature/orchid-literature.jsonl.gz"
    try:
        log = subprocess.run(["git", "-C", str(repo), "log", "--format=%H", "--", rel], capture_output=True, text=True,
                             check=True).stdout.split()
    except (OSError, subprocess.CalledProcessError):
        return None
    for c in log:
        b = subprocess.run(["git", "-C", str(repo), "show", "%s:%s" % (c, rel)], capture_output=True).stdout
        if sha256_bytes(b) == sample["input"]["sha256"]:
            return c
    return None


def redraw_matches(sample, sample_bytes, commit, repo=DATA_REPO):
    """Re-draw from the population snapshot at `commit` and compare bytes (input paths as recorded)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        lit = td / "data" / "literature"
        lit.mkdir(parents=True)
        for name in ("orchid-literature.jsonl.gz", "papers.json", "index.json"):
            b = subprocess.run(["git", "-C", str(repo), "show", "%s:literature/%s" % (commit, name)], capture_output=True).stdout
            (lit / name).write_bytes(b)
        if sha256_file(lit / "papers.json") != sample["input"]["papers_sha256"]:
            return False
        s = build(lit / "orchid-literature.jsonl.gz", lit / "papers.json", sample["design"]["seed"], sample["design"]["n"])
        return canonical(s).encode("utf-8") == sample_bytes


def cmd_check(a):
    bad, d = verify(a.out, a.input, a.papers)
    ch = [x for x in d if x["status"] == "changed"]
    if a.drift_json:
        a.drift_json.write_text(json.dumps({"input_sha256": sha256_file(a.input), "items": d}, sort_keys=True, indent=2,
                                           ensure_ascii=False) + "\n", encoding="utf-8")
    if bad:
        print("precision sample: " + "; ".join(bad[:10]) + (" ... (%d)" % len(bad) if len(bad) > 10 else ""), file=sys.stderr)
        return 1
    fz = json.loads((a.out / "frozen.json").read_text(encoding="utf-8"))
    mg = [x for x in d if x.get("merged_into")]
    print("precision sample frozen and intact (sample sha256 %s, %d items; population %s); all ids released (%d merged "
          "into a duplicate); "
          "%d item(s) with a released value changed since sampling%s" % (
              fz["sample_sha256"][:16], len(fz["items"]), fz["population"]["sha256"][:16], len(mg), len(ch),
              (": " + ", ".join("n%d %s [%s]" % (x["n"], x["id"], ",".join(x["questions"])) for x in ch[:8])
               + (" ..." if len(ch) > 8 else "")) if ch else ""))
    return 0


def cmd_freeze(a):
    sb = (a.out / "sample.json").read_bytes()
    sample = json.loads(sb)
    page, _ = render(sb.decode("utf-8"))
    if (a.out / "index.html").read_text(encoding="utf-8") != page:
        print("index.html does not embed sample.json; refusing to freeze", file=sys.stderr)
        return 1
    commit = population_commit(sample)
    reproduced = redraw_matches(sample, sb, commit) if (a.verify_draw and commit) else None
    if a.verify_draw and not reproduced:
        print("re-draw from the population snapshot (%s) does not reproduce sample.json" % (commit or "commit not found"),
              file=sys.stderr)
        return 1
    rec = freeze_record(a.out, sb, commit, reproduced)
    text = json.dumps(rec, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    fp = a.out / "frozen.json"
    if fp.exists():
        old = json.loads(fp.read_text(encoding="utf-8"))
        if old["files_sha256"] != rec["files_sha256"] or old["items"] != rec["items"]:
            print("frozen.json exists and differs: the sample is frozen; refusing to overwrite", file=sys.stderr)
            return 1
        if old.get("population", {}).get("draw_reproduced") and not reproduced:
            print("frozen.json is current")
            return 0
    fp.write_text(text, encoding="utf-8")
    print("frozen.json written: sample sha256 %s, %d items, population %s (orchid-data %s)%s" % (
        rec["sample_sha256"][:16], len(rec["items"]), rec["population"]["sha256"][:16], (commit or "?")[:9],
        ", re-draw reproduces the sample" if reproduced else ""))
    return 0


def cmd_draw(a):
    if a.out.resolve() == OUTDIR.resolve() or (a.out / "sample.json").exists() or (a.out / "frozen.json").exists():
        print("refusing to draw into %s: a sample exists there (samples are frozen; draw into a new directory)" % a.out,
              file=sys.stderr)
        return 1
    sample = build(a.input, a.papers, a.seed, a.n)
    text = canonical(sample)
    page, sha = render(text)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "sample.json").write_text(text, encoding="utf-8")
    (a.out / "index.html").write_text(page, encoding="utf-8")
    rec = freeze_record(a.out, text.encode("utf-8"))
    (a.out / "frozen.json").write_text(json.dumps(rec, sort_keys=True, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("new sample in %s: sha256 %s; %d items (frozen)" % (a.out, sha, len(sample["items"])))
    for s in sample["strata"]:
        print("  %-30s population %6d  sampled %3d" % (s["stratum"], s["population"], s["sampled"]))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("command", nargs="?", default="check", choices=["check", "freeze", "draw"])
    ap.add_argument("--input", type=Path, default=INPUT, help="literature JSONL.gz (default: %(default)s)")
    ap.add_argument("--papers", type=Path, default=PAPERS, help="papers.json (default: %(default)s)")
    ap.add_argument("--out", type=Path, default=OUTDIR, help="sample directory (default: %(default)s)")
    ap.add_argument("--seed", type=int, default=SEED, help="draw: seed")
    ap.add_argument("-n", type=int, default=N, help="draw: sample size (default: %(default)s)")
    ap.add_argument("--check", action="store_true", help="verify the frozen sample against the current export (default)")
    ap.add_argument("--drift-json", type=Path, help="check: write the per-item drift list (for the scorer) to this path")
    ap.add_argument("--verify-draw", action="store_true", help="freeze: require that a re-draw from the population "
                    "snapshot in the orchid-data git history reproduces sample.json byte for byte")
    a = ap.parse_args(argv)
    if a.check or a.command == "check":
        return cmd_check(a)
    return {"freeze": cmd_freeze, "draw": cmd_draw}[a.command](a)


if __name__ == "__main__":
    sys.exit(main())
