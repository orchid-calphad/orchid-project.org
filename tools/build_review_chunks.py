#!/usr/bin/env python3
"""Split the precision-sample review page into reviewer packets (deterministic).

Packet k gets every item whose sample position i satisfies i % K == k (so each packet mixes strata), plus the first
OVERLAP own items of packet (k+1) % K for inter-reviewer agreement. Each packet page embeds only its items but keeps the
full sample's SHA-256, so result files score against the one sample with tools/score_precision_sample.py.
Usage: build_review_chunks.py [--k 6] [--overlap 10] [--out DIR]"""
import argparse, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "review", "precision-sample")
ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--k", type=int, default=6); ap.add_argument("--overlap", type=int, default=10)
ap.add_argument("--out", default=os.path.join(SRC, "packets"))
a = ap.parse_args()
page = open(os.path.join(SRC, "index.html"), encoding="utf-8").read()
m = re.search(r'(<script type="application/json" id="sample-data">)(.*?)(</script>)', page, re.S)
sample = json.loads(m.group(2))
items = sample["items"]
own = [[it for i, it in enumerate(items) if i % a.k == k] for k in range(a.k)]
os.makedirs(a.out, exist_ok=True)
plan = []
for k in range(a.k):
    extra = own[(k + 1) % a.k][:a.overlap]
    sub = dict(sample, items=own[k] + extra,
               packet={"index": k + 1, "of": a.k, "own": len(own[k]), "overlap_with_packet": (k + 1) % a.k + 1, "overlap": len(extra)})
    body = json.dumps(sub, ensure_ascii=False, indent=1).replace("</", "<\\/")
    html = page[:m.start(2)] + body + page[m.end(2):]
    html = html.replace("<title>Precision sample review</title>", f"<title>ORCHID review packet {k+1} of {a.k}</title>", 1)
    name = f"ORCHID_review_packet_{k+1}.html"
    open(os.path.join(a.out, name), "w", encoding="utf-8").write(html)
    plan.append({"packet": k + 1, "file": name, "items": len(sub["items"]), "own": len(own[k]), "overlap_with_packet": (k + 1) % a.k + 1})
json.dump({"k": a.k, "overlap": a.overlap, "packets": plan}, open(os.path.join(a.out, "packets.json"), "w"), indent=1)
for p in plan: print(p)
