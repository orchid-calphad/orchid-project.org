#!/usr/bin/env python3
"""Unit tests of tools/quote_clip.py (export-time clipping of quotes longer than 500 characters).
usage: python3 tools/tests/test_quote_clip.py"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import quote_clip as QC  # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (("  " + str(detail)) if detail and not cond else ""))
    if not cond:
        fails.append(name)


filler = "The sample was annealed and then examined by optical microscopy and X-ray diffraction. "
short = "The eutectic is at 1273 K."
check("short quote untouched", QC.clip(short, ["1273"]) == (short, False, True))

for pos in ("start", "middle", "end"):
    parts = [filler * 4, "The eutectic temperature was found to be 1273 K in this work. ", filler * 4]
    if pos == "start":
        parts = [parts[1], filler * 8]
    elif pos == "end":
        parts = [filler * 8, parts[1]]
    q = "".join(parts).strip()
    assert len(q) > 500
    out, was, kept = QC.clip(q, ["1273"])
    check("value %s: clipped to <= 500 chars" % pos, was and len(out) <= 500, len(out))
    check("value %s: contiguous substring, value kept, no ellipsis" % pos, out in q and "1273" in out and kept
          and "..." not in out and "…" not in out)
    check("value %s: word boundaries" % pos, (out == q[:len(out)] or q[q.index(out) - 1].isspace())
          and (q.endswith(out) or q[q.index(out) + len(out)].isspace()))

q = ("Alpha. " * 40) + "A value of 12.5 kJ/mol was measured at 900 K. " + ("Beta. " * 60)
r = {"quote": q, "value_printed": "12.5", "temperature_printed": "900 K", "quote_page": 7, "flags": {}}
check("record clipped, flag, page kept, original length", QC.clip_record(r) and r["flags"]["quote_clipped"]
      and r["flags"]["quote_chars_original"] == len(q) and r["quote_page"] == 7 and "12.5" in r["quote"]
      and "900" in r["quote"] and len(r["quote"]) <= 500)
check("record under the limit: no flag", not QC.clip_record({"quote": short, "flags": {}}))
it = {"quote": filler * 12 + "Heated at 1073 K for 24 h. " + filler * 12, "temperature_printed": "1073 K"}
check("context item clipped around its printed temperature", QC.clip_context_item(it) and "1073" in it["quote"]
      and it["quote_clipped"] and len(it["quote"]) <= 500)
check("number match is whole-token (12 is not in 125)", QC._spans("a 125 b", "12") == [])
print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
