#!/usr/bin/env python3
"""Score the supplementary precision sample (result files downloaded from index.html or the packets).

Per sub-sample:
  S-REG, S-ESPEI  precision per question and for "all questions correct" (no question answered "no" and the
                  sub-sample's key question answered "yes": names for S-REG, value for S-ESPEI), per stratum cell, per
                  primary and secondary stratum, overall (plain) and population-weighted, with Wilson 95 % intervals.
                  "can't tell" and "n/a" are excluded from the denominator and counted.
  S-ASM           blind naming: the reviewer's name is compared with the hidden machine name (key.json) as a set of
                  phases after normalisation (Greek variants, +x subscripts, parentheses around terminal solutions,
                  L1 -> L, gas species -> G), then with the PI's notation equivalences (occurrence numbering of one phase,
                  alpha/beta prefixes printed without a transition line), exactly as the earlier blind sample's
                  scorer (open-tdb build/okamoto-pilot/tools/blind_score.py, copied below). Verdicts: agree,
                  agree_notation, disagree; cannot_tell / tracing_fault / source_undetermined are counted, not errors.
                  Precision = (agree + agree_notation) / named, also raw (agree / named).

Population weighting: strata are allocated with minimums, so the plain proportion over-represents small strata; the
weighted estimate uses stratum-cell population shares and its interval is a Wilson interval at the Kish effective n
(approximate, stated as such).

Which answer counts: every item has one owner packet (the kit's partition: partition.json named in manifest.json,
or the draw's packet plan in key.json when the kit was never re-partitioned; results from superseded packets are
refused). Precision uses the owner's answer (another
reviewer's answer only when the owner gave none); items answered by two reviewers (the ~10 % overlap) give the
inter-reviewer agreement: raw agreement and Cohen's kappa per question on items both judged yes/no, for S-ASM the
agreement of the two names (phases_equiv) and of the decision class, and per overlap pair of the partition (packet p
holds items owned by packet p+1) the share of matching judgements. Any number of packets.

Scores are computed against the frozen sample: sample.json and key.json must match manifest.json.

  score_supplement.py [--kit DIR] [--json OUT] RESULT.json [RESULT.json ...]
"""
import argparse
import collections
import hashlib
import json
import math
import re
import sys
from pathlib import Path

KIT = Path(__file__).resolve().parent.parent
Z = 1.959963984540054
QUESTIONS = {
    "S-REG": ("names", "structure", "tdb_definition", "example_dois"),
    "S-ESPEI": ("value", "unit", "conditions", "phases", "reference_state", "conversion", "uncertainty"),
}
KEY_QUESTION = {"S-REG": "names", "S-ESPEI": "value"}
NON_NAMES = ("cannot_tell", "tracing_fault", "source_undetermined")

# ---- phase-name normalisation: copied from open-tdb build/okamoto-pilot/tools/pd_fields.py (greek_norm) and
# blind_score.py (phases, phases_equiv), 2026-10-24, so that S-ASM is scored exactly as the earlier blind sample.
GREEK_VARIANTS = str.maketrans({"ϕ": "φ", "ϑ": "θ", "ϵ": "ε", "ϰ": "κ", "ϱ": "ρ", "ϖ": "π", "µ": "μ"})
PLUS_X = re.compile(r"(\d)\+([xy])(?![A-Za-z])|(\d)\+([xy])(?=[A-Z])")


def greek_norm(t):
    return PLUS_X.sub(lambda m: (m.group(1) or m.group(3)) + "₊" + (m.group(2) or m.group(4)), t.translate(GREEK_VARIANTS))


def phases(name):
    out = set()
    for t in greek_norm(name or "").replace(" ", "").split("+"):
        t = re.sub(r"^\((.*)\)$", r"\1", t)
        t = {"L1": "L", "Gas": "G", "gas": "G"}.get(t, t)
        if re.fullmatch(r"(O2|N2|H2|F2|Cl2)(\(G\))?", t):
            t = "G"
        if t:
            out.add(t)
    return frozenset(out)


GREEK = "αβγδεζηθικλμνξοπρστυφχψω"


def phases_equiv(name):
    out = collections.Counter()
    raw = (name or "").replace("alpha/", "α/").replace("beta", "β").replace(" ", "")
    for t in raw.split("+"):
        t = re.sub(r"^[%s]/[%s]" % (GREEK, GREEK), "", t)
        t = re.sub(r"#\d+$", "", t)
        t = re.sub(r"^(\(.*\))\d+$", r"\1", t)
        for x in phases(t):
            out[x] += 1
    return out


def asm_verdict(decision, machine):
    if decision in (None, "") or decision in NON_NAMES:
        return decision or "unanswered"
    if phases(decision) == phases(machine):
        return "agree"
    if phases_equiv(decision) == phases_equiv(machine):
        return "agree_notation"
    return "disagree"


# ---- statistics

def wilson(k, n, z=Z):
    if n <= 0:
        return (None, None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, c - hw), min(1.0, c + hw))


def est(k, n):
    p, lo, hi = wilson(k, n)
    return {"correct": k, "judged": n, "precision": p, "ci95": [lo, hi] if p is not None else None}


def weighted(per, pop):
    """per {cell: [k, n]} -> population-weighted precision over cells with n > 0, Wilson at Kish n_eff."""
    judged = {c: kn for c, kn in per.items() if kn[1] > 0}
    W = sum(pop.get(c, 0) for c in judged)
    if not W:
        return None
    wp = sum(pop[c] / W * kn[0] / kn[1] for c, kn in judged.items())
    neff = 1 / sum((pop[c] / W) ** 2 / kn[1] for c, kn in judged.items())
    _, lo, hi = wilson(wp * neff, neff)
    return {"precision": wp, "ci95": [lo, hi], "n_eff": neff, "population_covered": W}


def fmt(e):
    if e is None or e.get("precision") is None:
        return "%4d/%-4d      —" % ((e or {}).get("correct", 0), (e or {}).get("judged", 0))
    return "%4d/%-4d %6.1f%%  [%5.1f, %5.1f]" % (e["correct"], e["judged"], 100 * e["precision"], 100 * e["ci95"][0],
                                                 100 * e["ci95"][1])


def fmtw(w):
    if not w:
        return "—"
    return "%6.1f%%  [%5.1f, %5.1f]  n_eff %.1f" % (100 * w["precision"], 100 * w["ci95"][0], 100 * w["ci95"][1], w["n_eff"])


# ---- loading

def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_kit(kit):
    """-> (manifest, key, partition). The partition (owners, overlap pairs) comes from manifest.json /
    partition.json when the kit was re-partitioned, else from the packet plan of the draw (key.json)."""
    man = json.loads((kit / "manifest.json").read_text(encoding="utf-8"))
    for name, field in (("sample.json", "sample_sha256"), ("key.json", "key_sha256")):
        if sha256(kit / name) != man[field]:
            raise SystemExit("%s does not match manifest.json (the sample is frozen)" % (kit / name))
    key = json.loads((kit / "key.json").read_text(encoding="utf-8"))
    if man.get("partition"):
        pf = kit / man["partition"]["file"]
        if sha256(pf) != man["files_sha256"].get(man["partition"]["file"]):
            raise SystemExit("%s does not match manifest.json" % pf)
        part = json.loads(pf.read_text(encoding="utf-8"))
    else:
        part = {"id": None, "k": len(key["packets"]), "owner": {it["item"]: it["owner_packet"] for it in key["items"]},
                "packets": key["packets"]}
    return man, key, part


def load_results(paths, sample_sha, items, part=None):
    res = []
    for p in paths:
        r = json.loads(Path(p).read_text(encoding="utf-8"))
        if r.get("sample_sha256") != sample_sha:
            raise SystemExit("%s: sample_sha256 %s does not match the sample (%s)" % (p, r.get("sample_sha256"), sample_sha))
        if part and r.get("packet") is not None and r.get("partition") != part.get("id"):
            raise SystemExit("%s: packet %s of partition %s, but the kit's partition is %s (superseded packets are not "
                             "scored)" % (p, r.get("packet"), r.get("partition"), part.get("id")))
        if part and r.get("packet") is not None:
            allowed = {i for pk in part["packets"] if pk["packet"] == r["packet"] for i in pk["own"] + pk["overlap"]}
            stray = [a["item"] for a in r.get("answers") or [] if a.get("item") not in allowed]
            if stray:
                raise SystemExit("%s: items %s are not in packet %s" % (p, ", ".join(stray[:5]), r["packet"]))
        for a in r.get("answers") or []:
            it = items.get(a.get("item"))
            if it is None or it["kind"] != a.get("kind"):
                raise SystemExit("%s: answer %s is not in the sample or has the wrong kind" % (p, a.get("item")))
        r["_file"] = str(p)
        res.append(r)
    return res


def resolve(results, items, owner):
    """{item: [(result, answer)]} and the scored answer per item (owner packet first, then file order)."""
    by = collections.defaultdict(list)
    for r in results:
        for a in r["answers"]:
            by[a["item"]].append((r, a))
    scored = {}
    for i, lst in by.items():
        own = [ra for ra in lst if ra[0].get("packet") in (None, owner[i])]
        scored[i] = (own or lst)[0]
    return by, scored


# ---- scoring

def all_correct(a, kind):
    vals = [a.get(q) for q in QUESTIONS[kind]]
    if "no" in vals:
        return "no"
    if a.get(KEY_QUESTION[kind]) == "yes" and all(v in ("yes", "na") for v in vals):
        return "yes"
    return None


def score_questions(kind, answers, strata):
    """answers [(stratum, answer)] -> per question: overall, weighted, per cell, per primary, per secondary."""
    pop = {s["stratum"]: s["population"] for s in strata}
    out = {}
    for q in QUESTIONS[kind] + ("all",):
        per = collections.defaultdict(lambda: [0, 0])
        excluded = collections.Counter()
        for st, a in answers:
            v = all_correct(a, kind) if q == "all" else (a.get(q) if a.get(q) in ("yes", "no") else None)
            if v is None:
                excluded["incomplete" if q == "all" else (a.get(q) or "unanswered")] += 1
                continue
            per[st][1] += 1
            per[st][0] += v == "yes"
        tot = [sum(kn[0] for kn in per.values()), sum(kn[1] for kn in per.values())]
        marg = {}
        for lvl in (0, 1):
            m = collections.defaultdict(lambda: [0, 0])
            for st, kn in per.items():
                m[st.split("|")[lvl]][0] += kn[0]
                m[st.split("|")[lvl]][1] += kn[1]
            marg[("primary", "secondary")[lvl]] = {k: est(*kn) for k, kn in sorted(m.items())}
        out[q] = {"overall": est(*tot), "weighted": weighted(dict(per), pop), "excluded": dict(excluded),
                  "per_stratum": {s["stratum"]: est(*per.get(s["stratum"], [0, 0])) for s in strata},
                  "per_primary": marg["primary"], "per_secondary": marg["secondary"]}
    return out


def score_asm(rows, strata):
    pop = {s["stratum"]: s["population"] for s in strata}
    counts = collections.defaultdict(collections.Counter)
    for st, v in rows:
        counts[st][v] += 1
    res = {}
    for mode, ok in (("with_notation_equivalences", ("agree", "agree_notation")), ("raw", ("agree",))):
        per = {st: [sum(c[x] for x in ok), c["agree"] + c["agree_notation"] + c["disagree"]] for st, c in counts.items()}
        tot = [sum(kn[0] for kn in per.values()), sum(kn[1] for kn in per.values())]
        marg = {}
        for lvl in (0, 1):
            m = collections.defaultdict(lambda: [0, 0])
            for st, kn in per.items():
                m[st.split("|")[lvl]][0] += kn[0]
                m[st.split("|")[lvl]][1] += kn[1]
            marg[lvl] = {k: est(*kn) for k, kn in sorted(m.items())}
        res[mode] = {"overall": est(*tot), "weighted": weighted(per, pop),
                     "per_stratum": {s["stratum"]: est(*per.get(s["stratum"], [0, 0])) for s in strata},
                     "per_basis": marg[0], "per_boundary": marg[1]}
    res["verdicts"] = dict(sum(counts.values(), collections.Counter()))
    res["verdicts_per_stratum"] = {k: dict(v) for k, v in sorted(counts.items())}
    return res


def kappa(pairs):
    n = len(pairs)
    if not n:
        return None, None
    po = sum(x == y for x, y in pairs) / n
    cats = sorted({x for p in pairs for x in p})
    pe = sum((sum(x == c for x, _ in pairs) / n) * (sum(y == c for _, y in pairs) / n) for c in cats)
    return po, ((po - pe) / (1 - pe) if pe < 1 else None)


def agreement(by, items):
    out = {}
    multi = {i: lst for i, lst in by.items() if len({id(r) for r, _ in lst}) >= 2}
    for kind in ("S-REG", "S-ESPEI"):
        qa = {}
        for q in QUESTIONS[kind]:
            pairs, dis = [], []
            for i, lst in sorted(multi.items()):
                if items[i]["kind"] != kind:
                    continue
                x, y = lst[0][1].get(q), lst[1][1].get(q)
                if x in ("yes", "no") and y in ("yes", "no"):
                    pairs.append((x, y))
                    if x != y:
                        dis.append(i)
            po, k = kappa(pairs)
            qa[q] = {"n": len(pairs), "agreement": po, "kappa": k, "disagreements": dis}
        out[kind] = qa
    pairs_name, pairs_class, dis = [], [], []
    for i, lst in sorted(multi.items()):
        if items[i]["kind"] != "S-ASM":
            continue
        x, y = lst[0][1].get("decision"), lst[1][1].get("decision")
        if not x or not y:
            continue
        cx = x if x in NON_NAMES else "name"
        cy = y if y in NON_NAMES else "name"
        pairs_class.append((cx, cy))
        if cx == cy == "name":
            same = phases_equiv(x) == phases_equiv(y)
            pairs_name.append(same)
            if not same:
                dis.append(i)
    out["S-ASM"] = {"n_named_by_both": len(pairs_name),
                    "name_agreement": (sum(pairs_name) / len(pairs_name)) if pairs_name else None,
                    "decision_class": dict(zip(("agreement", "kappa"), kappa(pairs_class)), n=len(pairs_class)),
                    "disagreements": dis}
    out["n_items_answered_by_two_or_more"] = len(multi)
    return out


def pair_agreement(by, items, part):
    """Per overlap pair of the partition (packet p holds items owned by packet p+1): items answered in both packets
    and the share of matching judgements (yes/no questions; S-ASM names by phases_equiv)."""
    out = []
    for pk in part["packets"]:
        a_, b_ = pk["packet"], pk.get("overlap_from_packet")
        both, same, tot = 0, 0, 0
        for i in pk["overlap"]:
            x = next((ans for r, ans in by.get(i, []) if r.get("packet") == a_), None)
            y = next((ans for r, ans in by.get(i, []) if r.get("packet") == b_), None)
            if not (x and y):
                continue
            both += 1
            if items[i]["kind"] == "S-ASM":
                dx, dy = x.get("decision"), y.get("decision")
                if dx and dy:
                    tot += 1
                    same += (phases_equiv(dx) == phases_equiv(dy)) if dx not in NON_NAMES and dy not in NON_NAMES else dx == dy
            else:
                for q in QUESTIONS[items[i]["kind"]]:
                    if x.get(q) in ("yes", "no") and y.get(q) in ("yes", "no"):
                        tot += 1
                        same += x.get(q) == y.get(q)
        out.append({"packet": a_, "overlap_from_packet": b_, "overlap_items": len(pk["overlap"]), "answered_by_both": both,
                    "judgements": tot, "agreement": same / tot if tot else None})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results", nargs="+", type=Path)
    ap.add_argument("--kit", type=Path, default=KIT, help="kit directory with manifest.json, sample.json, key.json")
    ap.add_argument("--json", type=Path, help="also write the scores as JSON")
    a = ap.parse_args(argv)
    man, key, part = load_kit(a.kit)
    items = {it["item"]: it for it in key["items"]}
    results = load_results(a.results, man["sample_sha256"], items, part)
    by, scored = resolve(results, items, part["owner"])
    report = {"sample_sha256": man["sample_sha256"], "seed": key["design"]["seed"],
              "reviewers": [{"file": r["_file"], "reviewer": r.get("reviewer"), "packet": r.get("packet"),
                             "answers": len(r["answers"])} for r in results]}
    out = ["Supplementary sample %s (seed %s, %d items); %d result file(s)" % (
        man["sample_sha256"][:16], key["design"]["seed"], len(items), len(results))]
    for kind in ("S-REG", "S-ESPEI"):
        ans = [(items[i]["stratum"], a_) for i, (_, a_) in scored.items() if items[i]["kind"] == kind]
        sc = score_questions(kind, ans, key["strata"][kind])
        report[kind] = sc
        out += ["", "== %s: %d items answered ==" % (kind, len(ans)),
                "  %-16s %-34s %s" % ("question", "correct/judged  precision [95% CI]", "population-weighted")]
        for q in QUESTIONS[kind] + ("all",):
            ex = ", ".join("%s %d" % kv for kv in sorted(sc[q]["excluded"].items()))
            out.append("  %-16s %-34s %s%s" % ("all questions" if q == "all" else q, fmt(sc[q]["overall"]),
                                               fmtw(sc[q]["weighted"]), ("   excluded: " + ex) if ex else ""))
        out.append("  per stratum (all questions):")
        for s in key["strata"][kind]:
            out.append("    %-48s %s" % (s["stratum"], fmt(sc["all"]["per_stratum"][s["stratum"]])))
        for lvl in ("per_primary", "per_secondary"):
            out.append("  %s (all questions): %s" % (lvl.replace("_", " "), "; ".join(
                "%s %s" % (k, fmt(v).strip()) for k, v in sc["all"][lvl].items())))
    rows = []
    for i, (_, a_) in scored.items():
        it = items[i]
        if it["kind"] == "S-ASM":
            rows.append((it["stratum"], asm_verdict(a_.get("decision"), it["secret"]["machine_name"])))
    sa = score_asm(rows, key["strata"]["S-ASM"])
    sa["items"] = sorted(
        [{"item": i, "stratum": items[i]["stratum"], "machine": items[i]["secret"]["machine_name"],
          "reviewer": a_.get("decision"), "confidence": a_.get("confidence"), "note": a_.get("note"),
          "verdict": asm_verdict(a_.get("decision"), items[i]["secret"]["machine_name"])}
         for i, (_, a_) in scored.items() if items[i]["kind"] == "S-ASM"], key=lambda x: x["item"])
    report["S-ASM"] = sa
    out += ["", "== S-ASM (blind naming): %d items answered; verdicts %s ==" % (len(rows), json.dumps(sa["verdicts"], sort_keys=True))]
    for mode in ("with_notation_equivalences", "raw"):
        m = sa[mode]
        out.append("  %-28s %-34s %s" % (mode, fmt(m["overall"]), fmtw(m["weighted"])))
    out.append("  per stratum (with notation equivalences):")
    for s in key["strata"]["S-ASM"]:
        out.append("    %-48s %s" % (s["stratum"], fmt(sa["with_notation_equivalences"]["per_stratum"][s["stratum"]])))
    for lvl in ("per_basis", "per_boundary"):
        out.append("  %s: %s" % (lvl.replace("_", " "), "; ".join(
            "%s %s" % (k, fmt(v).strip()) for k, v in sa["with_notation_equivalences"][lvl].items())))
    for x in sa["items"]:
        if x["verdict"] in ("disagree", "agree_notation"):
            out.append("  %s %s machine: %s | reviewer: %s" % (x["verdict"].upper(), x["item"], x["machine"], x["reviewer"]))
    ag = agreement(by, items)
    ag["pairs"] = pair_agreement(by, items, part)
    report["agreement"] = ag
    report["partition"] = part.get("id")
    out += ["", "== Inter-reviewer agreement: %d items answered by two or more reviewers ==" % ag["n_items_answered_by_two_or_more"]]
    for kind in ("S-REG", "S-ESPEI"):
        for q, v in ag[kind].items():
            out.append("  %-8s %-16s n %3d  agreement %s  kappa %s  disagreements: %s" % (
                kind, q, v["n"], "  —  " if v["agreement"] is None else "%5.1f%%" % (100 * v["agreement"]),
                "  — " if v["kappa"] is None else "%5.2f" % v["kappa"], ", ".join(v["disagreements"]) or "none"))
    s = ag["S-ASM"]
    out.append("  S-ASM    names (both named) n %d  agreement %s; decision class n %d agreement %s kappa %s; disagreements: %s" % (
        s["n_named_by_both"], "—" if s["name_agreement"] is None else "%.1f%%" % (100 * s["name_agreement"]),
        s["decision_class"]["n"], "—" if s["decision_class"]["agreement"] is None else "%.1f%%" % (100 * s["decision_class"]["agreement"]),
        "—" if s["decision_class"]["kappa"] is None else "%.2f" % s["decision_class"]["kappa"], ", ".join(s["disagreements"]) or "none"))
    for pr in ag["pairs"]:
        out.append("  packets %d <- %d: %d overlap items, %d answered by both; agreement over yes/no questions and names %s" % (
            pr["packet"], pr["overlap_from_packet"], pr["overlap_items"], pr["answered_by_both"],
            "—" if pr["agreement"] is None else "%.1f%% (%d judgements)" % (100 * pr["agreement"], pr["judgements"])))
    print("\n".join(out))
    if a.json:
        a.json.write_text(json.dumps(report, sort_keys=True, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
