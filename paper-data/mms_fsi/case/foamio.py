"""Minimal ASCII OpenFOAM field reader and exact-solution evaluation helpers."""
import os
import re
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def _parse_list(txt, start):
    """Parse 'N ( ... )' starting at index start; return (array, end)."""
    m = re.compile(r"\s*(\d+)\s*\(").match(txt, start)
    n = int(m.group(1))
    i = m.end()
    depth = 1
    j = i
    while depth:
        c = txt[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        j += 1
    body = txt[i:j - 1]
    if "(" in body:
        vals = np.array([[float(v) for v in t.split()] for t in re.findall(r"\(([^()]*)\)", body)])
    else:
        vals = np.array([float(v) for v in body.split()])
    assert len(vals) == n, (len(vals), n)
    return vals, j


def read_field(path, patch=None):
    """Internal field (patch=None) or a patch's value list."""
    txt = open(path).read()
    if patch is None:
        k = txt.index("internalField")
        seg = txt[k:]
    else:
        b = txt.index("boundaryField")
        k = re.compile(r"\b" + re.escape(patch) + r"\b\s*\{").search(txt, b).end()
        seg = txt[k:]
        key = "value" if re.search(r"\bvalue\b", seg[:seg.index("}")] if "}" in seg else seg) else "gradient"
        seg = seg[re.search(r"\b" + key + r"\b", seg).start():]
    m = re.match(r"(internalField|value|gradient)\s+(uniform|nonuniform)\s+(List<\w+>)?", seg)
    if m.group(2) == "uniform":
        v = re.match(r"(internalField|value|gradient)\s+uniform\s+(\([^)]*\)|[-\d.eE+]+)", seg).group(2)
        return np.array([float(q) for q in v.strip("()").split()])
    vals, _ = _parse_list(seg, m.end())
    return vals


def exact(kind, pts, exe="evalExact"):
    """Evaluate the exact solution with ../evalExact for rows of points."""
    inp = "\n".join(f"{kind} " + " ".join(f"{q:.17g}" for q in p) for p in pts) + "\n"
    out = subprocess.run([os.path.join(HERE, exe)], input=inp, capture_output=True, text=True).stdout
    return np.array([[float(v) for v in l.split()] for l in out.strip().splitlines()])
