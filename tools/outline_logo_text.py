#!/usr/bin/env python3
"""Convert the <text> elements of the ORCHID logo SVGs to outlined <path>s using the IBM Plex font files in
assets/logo/fonts (SIL OFL). Advance widths only (no kerning), plus the letter-spacing given in LOGO_SPEC.md.
Writes *-outlined.svg next to each source; the live-text sources are kept."""
import re, sys, os
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
D = os.path.join(os.path.dirname(__file__), "..", "assets", "logo")
FONTS = {"sans600": "IBMPlexSans-SemiBold.ttf", "sans400": "IBMPlexSans-Regular.ttf", "serifi": "IBMPlexSerif-Italic.ttf"}
cache = {}
def path_for(fontkey, size, x, y, text, tracking=0.0):
    f = cache.setdefault(fontkey, TTFont(os.path.join(D, "fonts", FONTS[fontkey])))
    gs, cmap, upm = f.getGlyphSet(), f.getBestCmap(), f["head"].unitsPerEm
    s = size / upm; out = []; cx = x
    for ch in text:
        g = cmap.get(ord(ch))
        if g is None: continue
        pen = SVGPathPen(gs)
        gs[g].draw(TransformPen(pen, (s, 0, 0, -s, cx, y)))
        d = pen.getCommands()
        if d: out.append(d)
        cx += gs[g].width * s + tracking
    return " ".join(out), cx - x
def fmt(d):  # round numbers to 2 decimals
    return re.sub(r"-?\d+\.\d+", lambda m: f"{float(m.group()):.2f}".rstrip("0").rstrip("."), d)
def outline(src, specs):
    s = open(os.path.join(D, src)).read()
    for pat, (key, size, x, y, txt, tr, fill) in specs.items():
        d, w = path_for(key, size, x, y, txt, tr)
        s = s.replace(pat, f'<path d="{fmt(d)}" fill="{fill}"/>')
    s = re.sub(r"\s*<style>.*?</style>", "", s, flags=re.S)
    s = s.replace('<g font-family="\'IBM Plex Serif\', Georgia, serif" font-style="italic" font-size="12" fill="#1D1B20">', "<g>")
    out = src.replace(".svg", "-outlined.svg"); open(os.path.join(D, out), "w").write(s); print(out, "written")
lab = {'<text x="89" y="7">A</text>': ("serifi", 12, 89, 7, "A", 0, "#1D1B20"),
       '<text x="0" y="168">B</text>': ("serifi", 12, 0, 168, "B", 0, "#1D1B20"),
       '<text x="176" y="168">C</text>': ("serifi", 12, 176, 168, "C", 0, "#1D1B20")}
outline("orchid-mark.svg", {k.replace("    ", ""): v for k, v in lab.items()})
lock = dict(lab)
lock['<text class="wm" x="214" y="92">ORCHID</text>'] = ("sans600", 58, 214, 92, "ORCHID", 1.74, "#1D1B20")
lock['<text class="tg" x="214" y="120">Open Repository for CALPHAD Hypotheses,</text>'] = ("sans400", 14, 214, 120, "Open Repository for CALPHAD Hypotheses,", 0, "#5E5A66")
lock['<text class="tg" x="214" y="139">Intelligence, and Discovery</text>'] = ("sans400", 14, 214, 139, "Intelligence, and Discovery", 0, "#5E5A66")
outline("orchid-lockup.svg", lock)
