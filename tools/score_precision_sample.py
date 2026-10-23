#!/usr/bin/env python3
"""Score precision-sample review results (files downloaded from review/precision-sample/index.html).

For each result file, and for all files pooled: precision per field (value, unit, evidence_class, phases,
citation) and for "all fields correct", overall, per stratum and population-weighted, with Wilson 95 %
intervals. "can't tell" and "n/a" answers are excluded from the denominator and counted separately.
With two result files, inter-reviewer agreement per field: items both answered with yes/no, raw agreement
and Cohen's kappa.

Population weighting: the sample is stratified with a minimum per stratum, so the plain sample proportion
over-represents small strata. The weighted estimate uses stratum population shares; its interval is a Wilson
interval at the Kish effective sample size (an approximation, stated as such in the output).

Scores are always computed against the FROZEN sample (review/precision-sample/sample.json, verified against
frozen.json; see tools/build_precision_sample.py). The literature release moves on after sampling, so the report
also says, per answered item, whether the current release still carries the value the reviewer judged
(`unchanged`), carries a different value for one of the reviewed fields (`changed`, with the review questions
affected), or no longer releases the record (`missing`); `merged_into` names the surviving record when the
de-duplication merged the sampled one. A judgement on a `changed` item describes the sampled value, not the release.

Standard library only.
"""
import argparse
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_precision_sample as bps  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SAMPLE = ROOT / "review" / "precision-sample" / "sample.json"
FIELDS = ("value", "unit", "evidence_class", "phases", "citation")
Z = 1.959963984540054


def wilson(k, n, z=Z):
    if n <= 0:
        return (None, None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, c - hw), min(1.0, c + hw))


def fmt(k, n):
    p, lo, hi = wilson(k, n)
    if p is None:
        return "%4d/%-4d      —" % (k, n)
    return "%4d/%-4d %6.1f%%  [%5.1f, %5.1f]" % (k, n, 100 * p, 100 * lo, 100 * hi)


def load_results(path, sample_sha):
    r = json.loads(Path(path).read_text(encoding="utf-8"))
    if r.get("sample_sha256") != sample_sha:
        raise SystemExit("%s: sample_sha256 %s does not match the sample (%s)" % (path, r.get("sample_sha256"), sample_sha))
    return r


def binary(a, field):
    v = a.get(field)
    return v if v in ("yes", "no") else None


def all_correct(a):
    """yes when every answered field is yes and the value is yes; no when any field is no; else None."""
    vals = [a.get(f) for f in FIELDS]
    if "no" in vals:
        return "no"
    if a.get("value") == "yes" and all(v in ("yes", "na") for v in vals):
        return "yes"
    return None


def score(answers, strata_pop):
    rows = {}
    for field in FIELDS + ("all",):
        tot = [0, 0]
        per = defaultdict(lambda: [0, 0])
        excluded = defaultdict(int)
        for a in answers:
            v = all_correct(a) if field == "all" else binary(a, field)
            if v is None:
                excluded[(a.get(field) or "unanswered") if field != "all" else "incomplete"] += 1
                continue
            tot[1] += 1
            per[a["stratum"]][1] += 1
            if v == "yes":
                tot[0] += 1
                per[a["stratum"]][0] += 1
        # population-weighted estimate over strata with at least one judged item
        judged = {s: kn for s, kn in per.items() if kn[1] > 0}
        W = sum(strata_pop[s] for s in judged)
        wp = None
        neff = 0
        if W:
            wp = sum(strata_pop[s] / W * kn[0] / kn[1] for s, kn in judged.items())
            neff = 1 / sum((strata_pop[s] / W) ** 2 / kn[1] for s, kn in judged.items())
        rows[field] = {"total": tot, "per_stratum": dict(per), "excluded": dict(excluded), "weighted": wp, "n_eff": neff}
    return rows


def report(title, answers, strata, out):
    pop = {s["stratum"]: s["population"] for s in strata}
    rows = score(answers, pop)
    out.append("")
    out.append("== %s: %d answered items ==" % (title, len(answers)))
    out.append("Precision with Wilson 95% intervals (excluding can't tell / n/a)")
    out.append("  %-15s %s" % ("field", "correct/judged   precision  [95% CI]"))
    for f in FIELDS + ("all",):
        r = rows[f]
        ex = ", ".join("%s %d" % (k, v) for k, v in sorted(r["excluded"].items()))
        out.append("  %-15s %s%s" % ("all fields" if f == "all" else f, fmt(*r["total"]), ("   excluded: " + ex) if ex else ""))
    out.append("Population-weighted (stratum shares; interval: Wilson at Kish effective n, approximate)")
    for f in FIELDS + ("all",):
        r = rows[f]
        if r["weighted"] is None:
            continue
        neff = r["n_eff"]
        _, lo, hi = wilson(r["weighted"] * neff, neff)
        out.append("  %-15s %6.1f%%  [%5.1f, %5.1f]  n_eff %.1f" % ("all fields" if f == "all" else f, 100 * r["weighted"], 100 * lo, 100 * hi, neff))
    out.append("Per stratum (value / all fields)")
    for s in strata:
        k = s["stratum"]
        v = rows["value"]["per_stratum"].get(k, [0, 0])
        a = rows["all"]["per_stratum"].get(k, [0, 0])
        out.append("  %-30s value %s   all %s" % (k, fmt(*v), fmt(*a)))
    return rows


def agreement(r1, r2, out):
    a1 = {a["id"]: a for a in r1["answers"]}
    a2 = {a["id"]: a for a in r2["answers"]}
    common = sorted(set(a1) & set(a2))
    out.append("")
    out.append("== Inter-reviewer agreement: %s vs %s, %d items answered by both ==" % (r1.get("reviewer"), r2.get("reviewer"), len(common)))
    for f in FIELDS:
        pairs = [(a1[i].get(f), a2[i].get(f)) for i in common]
        pairs = [(x, y) for x, y in pairs if x in ("yes", "no") and y in ("yes", "no")]
        n = len(pairs)
        if not n:
            out.append("  %-15s no items judged yes/no by both" % f)
            continue
        po = sum(x == y for x, y in pairs) / n
        p1 = sum(x == "yes" for x, _ in pairs) / n
        p2 = sum(y == "yes" for _, y in pairs) / n
        pe = p1 * p2 + (1 - p1) * (1 - p2)
        kappa = (po - pe) / (1 - pe) if pe < 1 else float("nan")
        dis = [i for i in common if a1[i].get(f) in ("yes", "no") and a2[i].get(f) in ("yes", "no") and a1[i].get(f) != a2[i].get(f)]
        out.append("  %-15s n %3d  agreement %5.1f%%  kappa %s  disagreements: %s" % (
            f, n, 100 * po, "  n/a" if math.isnan(kappa) else "%5.2f" % kappa, ", ".join(dis[:12]) + (" ..." if len(dis) > 12 else "") or "none"))


def release_status(sample_path, input_path, papers_path):
    """{id: drift entry} for the frozen sample against the current release (verifies the freeze record first)."""
    fp = Path(sample_path).parent / "frozen.json"
    if not fp.exists():
        raise SystemExit("%s missing: the sample is not frozen (tools/build_precision_sample.py freeze)" % fp)
    frozen = json.loads(fp.read_text(encoding="utf-8"))
    sha = hashlib.sha256(Path(sample_path).read_bytes()).hexdigest()
    if sha != frozen["sample_sha256"]:
        raise SystemExit("%s (sha256 %s) is not the frozen sample (%s)" % (sample_path, sha[:16], frozen["sample_sha256"][:16]))
    recs, papers = bps.load_release(input_path, papers_path)
    return {x["id"]: x for x in bps.drift(frozen, recs, papers)}


def drift_report(results, status, out):
    answered = sorted({a["id"]: a["n"] for r in results for a in r["answers"]}.items(), key=lambda kv: kv[1])
    c = Counter(status[i]["status"] for i, _ in answered)
    out.append("")
    out.append("== Release check: does the current release still carry the reviewed value? ==")
    out.append("  %d answered items: %s" % (len(answered), ", ".join("%s %d" % kv for kv in sorted(c.items())) or "none"))
    for i, n in answered:
        x = status[i]
        line = "  n%-4d %-22s %s" % (n, i, x["status"])
        if x.get("questions"):
            line += " (%s)" % ", ".join(x["questions"])
        if x.get("merged_into"):
            line += " merged into %s" % x["merged_into"]
        out.append(line)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results", nargs="+", type=Path, help="result JSON file(s) downloaded from the review page")
    ap.add_argument("--sample", type=Path, default=SAMPLE, help="sample.json (default: %(default)s)")
    ap.add_argument("--json", type=Path, help="also write the scores as canonical JSON to this path")
    ap.add_argument("--input", type=Path, default=bps.INPUT, help="current literature export (default: %(default)s)")
    ap.add_argument("--papers", type=Path, default=bps.PAPERS, help="current papers.json (default: %(default)s)")
    a = ap.parse_args(argv)
    status = release_status(a.sample, a.input, a.papers)
    sample_bytes = a.sample.read_bytes()
    sha = hashlib.sha256(sample_bytes).hexdigest()
    sample = json.loads(sample_bytes)
    ids = {it["record"]["id"]: it["stratum"] for it in sample["items"]}
    results = [load_results(p, sha) for p in a.results]
    for p, r in zip(a.results, results):
        for ans in r["answers"]:
            if ids.get(ans["id"]) != ans.get("stratum"):
                raise SystemExit("%s: answer %s is not in the sample or has the wrong stratum" % (p, ans["id"]))
    out = ["Sample %s (%d items, seed %s)" % (sha[:16], len(ids), sample["design"]["seed"])]
    scores = {}
    for p, r in zip(a.results, results):
        scores[str(p)] = report("%s (reviewer: %s, %s)" % (p.name, r.get("reviewer"), r.get("timestamp")), r["answers"], sample["strata"], out)
    if len(results) > 1:
        pooled = [ans for r in results for ans in r["answers"]]
        scores["pooled"] = report("Pooled over %d reviewers (item-answers, not independent)" % len(results), pooled, sample["strata"], out)
    if len(results) == 2:
        agreement(results[0], results[1], out)
    drift_report(results, status, out)
    print("\n".join(out))
    if a.json:
        answered = sorted({ans["id"] for r in results for ans in r["answers"]})
        rel = {i: {k: v for k, v in status[i].items() if k in ("status", "questions", "merged_into")} for i in answered}
        a.json.write_text(json.dumps({"sample_sha256": sha, "scores": scores, "release_check": rel}, sort_keys=True,
                                     indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
