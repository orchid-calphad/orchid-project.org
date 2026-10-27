#!/usr/bin/env python3
"""Clip long quotes at EXPORT time (PI decision 2026-10-24: short quotes, never whole paragraphs).

A released `quote` is a short verbatim excerpt kept so that a value can be verified on its page. Quotes longer than
MAX_QUOTE_CHARS (500) are replaced by a contiguous substring of the original quote of at most 500 characters that
still contains the printed value:

* the window starts and ends on word boundaries; nothing is inserted (no ellipsis), nothing is reworded;
* it covers every number of the printed value (`value_printed`, else `temperature_printed`), then as many of the other
  printed numbers of the record as fit, and prefers to start at a sentence start and to end at a sentence end;
* the record gets `flags.quote_clipped = true`, `flags.quote_chars_original` and `flags.quote_original_sha256` (first 16 hex characters of the SHA-256 of the original quote: values read from one long sentence share it, so passages can still be counted); `quote_page` is untouched.

The source record store is never changed: the original quote stays where the reader put it, only the exported copy is
clipped. A quote whose printed value cannot be placed inside 500 characters is still clipped around the first
occurrence of the value and flagged `flags.quote_value_kept = false` (none in release 2026.10).

Standard library only.
"""
import hashlib
import re

MAX_QUOTE_CHARS = 500
_NUM = re.compile(r"(?<![\d.])\d+(?:\.\d+)?")
_SENT_END = re.compile(r"[.!?]$")


def numbers(*printed):
    """Distinct numeric tokens of printed strings (order kept)."""
    out = []
    for p in printed:
        if p is None or p == "":
            continue
        for t in _NUM.findall(str(p)):
            if t not in out:
                out.append(t)
    return out


def _spans(text, tok):
    return [(m.start(), m.end()) for m in re.finditer(r"(?<![\d.])" + re.escape(tok) + r"(?!\d)", text)]


def clip(text, required=(), optional=(), limit=MAX_QUOTE_CHARS):
    """-> (clipped_text, was_clipped, value_kept). `required` / `optional`: numeric tokens (see numbers())."""
    if text is None or len(text) <= limit:
        return text, False, True
    starts = [0] + [m.end() for m in re.finditer(r"\s+", text) if m.end() < len(text)]
    ends = [m.start() for m in re.finditer(r"\s+", text) if m.start() > 0] + [len(text)]
    req = [_spans(text, t) for t in required]
    req = [s for s in req if s]
    opt = [_spans(text, t) for t in optional]
    opt = [s for s in opt if s]
    best = None
    ei = 0
    for s in starts:
        # last word end with e - s <= limit
        e = None
        for x in ends:
            if x - s <= limit and x > s:
                e = x
        if e is None:
            continue
        cov_r = sum(1 for sp in req if any(a >= s and b <= e for a, b in sp))
        cov_o = sum(1 for sp in opt if any(a >= s and b <= e for a, b in sp))
        sent = (s == 0 or bool(re.search(r"[.!?]\s*$", text[:s]))) + bool(_SENT_END.search(text[s:e]))
        centre = 0.0
        if req:
            first = req[0][0]
            centre = -abs((s + e) / 2 - (first[0] + first[1]) / 2)
        score = (cov_r, cov_o, e - s >= 150, sent, e - s, centre)
        if best is None or score > best[0]:
            best = (score, s, e)
    if best is None:   # a single word longer than the limit: cannot happen for prose; keep the head
        return text[:limit], True, False
    _, s, e = best
    return text[s:e], True, len(req) == len([t for t in required]) and best[0][0] == len(req)


def clip_record(r):
    """Clip r['quote'] in place (running-text literature record); returns True if clipped."""
    q = r.get("quote")
    if not q or len(q) <= MAX_QUOTE_CHARS:
        return False
    vp = r.get("value_printed")
    req = numbers(vp) or numbers(r.get("temperature_printed")) or numbers(r.get("value") if vp in (None, "") else None)
    opt = [t for t in numbers(vp, r.get("temperature_printed")) if t not in req]
    new, was, kept = clip(q, req, opt)
    if not was:
        return False
    r["quote"] = new
    fl = r.setdefault("flags", {})
    fl["quote_clipped"] = True
    fl["quote_chars_original"] = len(q)
    fl["quote_original_sha256"] = hashlib.sha256(q.encode("utf-8")).hexdigest()[:16]   # identifies the passage (several values share one)
    if not kept:
        fl["quote_value_kept"] = False
    return True


def clip_context_item(it):
    """Same for an experimental-context item (no printed value; its printed temperature / duration anchor the window)."""
    q = it.get("quote")
    if not q or len(q) <= MAX_QUOTE_CHARS:
        return False
    req = numbers(it.get("temperature_printed")) or numbers(it.get("duration"))
    opt = [t for t in numbers(it.get("duration"), it.get("purity"), it.get("stated_uncertainty")) if t not in req]
    new, was, kept = clip(q, req, opt)
    if not was:
        return False
    it["quote"] = new
    it["quote_clipped"] = True
    it["quote_chars_original"] = len(q)
    return True
