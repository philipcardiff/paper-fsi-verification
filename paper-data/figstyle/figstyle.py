"""Helpers for the paper figures: load paper.mplstyle with Computer Modern text.

    import sys; sys.path.insert(0, "<repo>/paper-data/figstyle")
    import figstyle; figstyle.use()

paper.mplstyle asks for Computer Modern (CMU Serif). Where it is not installed,
use() takes the CMU Serif OpenType files of the local TeX installation
(kpsewhich cmunrm.otf ...), converts their outlines to TrueType once (cached in
~/.cache/paper-figstyle) and registers them with matplotlib, so the PDFs embed
clean TrueType subsets. Without TeX it falls back to the style's DejaVu Serif.
"""
import os
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STYLE = HERE / "paper.mplstyle"
TEXTWIDTH = 5.16          # in, full text width (372 pt)
HALFWIDTH = 2.5           # in
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
CACHE = Path(os.environ.get("PAPER_FIGSTYLE_CACHE", Path.home() / ".cache" / "paper-figstyle"))
CMU = {"cmunrm": "regular", "cmunti": "italic", "cmunbx": "bold", "cmunbi": "bold italic"}


def _otf_to_ttf(src, dst):
    from fontTools.ttLib import TTFont, newTable
    from fontTools.pens.cu2quPen import Cu2QuPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    f = TTFont(src)
    gs = f.getGlyphSet()
    glyf = newTable("glyf"); glyf.glyphOrder = f.getGlyphOrder(); glyf.glyphs = {}
    for name in f.getGlyphOrder():
        pen = TTGlyphPen(gs)
        gs[name].draw(Cu2QuPen(pen, 1.0, reverse_direction=True))
        glyf.glyphs[name] = pen.glyph()
    f["glyf"] = glyf
    f["loca"] = newTable("loca")
    f["maxp"].tableVersion = 0x00010000
    for k in ("maxZones", "maxTwilightPoints", "maxStorage", "maxFunctionDefs", "maxInstructionDefs",
              "maxStackElements", "maxSizeOfInstructions", "maxComponentElements"):
        setattr(f["maxp"], k, 0)
    f["maxp"].maxZones = 1
    f["maxp"].maxComponentDepth = 0
    f["head"].glyphDataFormat = 0
    post = f["post"]; post.formatType = 2.0; post.extraNames = []; post.mapping = {}
    del f["CFF "]
    if "VORG" in f:
        del f["VORG"]
    f.sfntVersion = "\x00\x01\x00\x00"
    f.save(dst)


def _register_cmu():
    if any(e.name == "CMU Serif" for e in fm.fontManager.ttflist):
        return True
    CACHE.mkdir(parents=True, exist_ok=True)
    for stem in CMU:
        ttf = CACHE / f"{stem}.ttf"
        if not ttf.exists():
            try:
                otf = subprocess.run(["kpsewhich", stem + ".otf"], capture_output=True,
                                     text=True).stdout.strip()
            except OSError:
                return False
            if not otf:
                return False
            _otf_to_ttf(otf, ttf)
        fm.fontManager.addfont(str(ttf))
    return True


def use():
    """Apply paper.mplstyle, with Computer Modern text where available."""
    import logging
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    _register_cmu()
    plt.style.use(str(STYLE))
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


def panel_label(ax, s, x=0.0, y=1.02, **kw):
    """Panel label '(a)' at the top left, above the axes."""
    kw.setdefault("ha", "left"); kw.setdefault("va", "bottom")
    return ax.text(x, y, s, transform=ax.transAxes, **kw)
