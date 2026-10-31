#!/usr/bin/env python3
"""Minimal FEBio .xplt (format 0x35, uncompressed) reader extracted unchanged from the project tool tools/cavity_xplt_audit.py (reader part only; the cavity-specific analysis is not included).
Chunk = u32 tag, u32 size, payload; format from FEBio v4.13 FEBioPlot/FEBioPlotFile.{h,cpp}.  read_xplt(path) -> dict(dic, states[{time, flag, vars}], ids, X0, domains)."""
import sys, os, json, struct, numpy as np
BR = {0x01000000, 0x01010000, 0x01020000, 0x01021000, 0x01023000, 0x01024000, 0x01025000, 0x01026000, 0x01020001, 0x01040000, 0x01041000, 0x01041100, 0x01042000, 0x01042100, 0x01042101, 0x01042200,
      0x01043000, 0x01043100, 0x01043200, 0x01044000, 0x01044100, 0x01045000, 0x01045100, 0x01046000, 0x01046100, 0x01047000, 0x01047100, 0x01048000, 0x01048100, 0x01048200, 0x01050000,
      0x02000000, 0x02010000, 0x02020000, 0x02020001, 0x02020003, 0x02020300, 0x02020400, 0x02020500, 0x02020600, 0x02030000, 0x02040000, 0x01051000, 0x01052000}
# tags that are branches in the writer; PLT_STATE_VAR_DATA (0x02020003) holds either one raw float block (node data: tag 0) or (domain id, floats) sub-chunks.
def chunks(b, i, end):
    while i + 8 <= end:
        tag, sz = struct.unpack_from("<II", b, i); yield tag, i + 8, sz; i += 8 + sz
def read_xplt(path):
    b = open(path, "rb").read(); assert b[:4] == b"BEF\0"
    out = {"dic": {"node": [], "domain": []}, "states": [], "domains": []}; trunc = False
    def dic(i, sz, key):
        for t, p, s in chunks(b, i, i + sz):
            it = {}
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01020002: it["type"] = struct.unpack_from("<I", b, p2)[0]
                elif t2 == 0x01020003: it["fmt"] = struct.unpack_from("<I", b, p2)[0]
                elif t2 == 0x01020004: it["name"] = b[p2:p2 + s2].split(b"\0")[0].decode()
            out["dic"][key].append(it)
    NV = {0: 1, 1: 3, 2: 6, 3: 3, 4: 21, 5: 9}
    for t, p, s in chunks(b, 4, len(b)):
        if t == 0x01000000:
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01020000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x01023000: dic(p3, s3, "node")
                        if t3 == 0x01024000: dic(p3, s3, "domain")
        elif t == 0x01040000:
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01041000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x01041200: a = np.frombuffer(b, np.float32, s3 // 4, p3).reshape(-1, 4); out["ids"] = a[:, 0].copy().view(np.int32); out["X0"] = a[:, 1:].astype(float)
                if t2 == 0x01042000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 != 0x01042100: continue
                        d = {"elems": []}
                        for t4, p4, s4 in chunks(b, p3, p3 + s3):
                            if t4 == 0x01042101:
                                for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                    if t5 == 0x01042102: d["etype"] = struct.unpack_from("<I", b, p5)[0]
                                    if t5 == 0x01042103: d["part"] = struct.unpack_from("<I", b, p5)[0]
                            if t4 == 0x01042200:
                                for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                    d["elems"].append(np.frombuffer(b, np.int32, s5 // 4, p5).copy())
                        out["domains"].append(d)
        elif t == 0x02000000:
            st = {"vars": {}, "complete": True}
            if p + s > len(b): trunc = True; break
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x02010000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x02010002: st["time"] = struct.unpack_from("<f", b, p3)[0]
                        if t3 == 0x02010003: st["flag"] = struct.unpack_from("<I", b, p3)[0]
                if t2 == 0x02020000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        kind = {0x02020300: "node", 0x02020400: "domain"}.get(t3)
                        if kind is None: continue
                        for t4, p4, s4 in chunks(b, p3, p3 + s3):   # PLT_STATE_VARIABLE
                            vid = None
                            for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                if t5 == 0x02020002: vid = struct.unpack_from("<I", b, p5)[0]
                                if t5 == 0x02020003:
                                    it = out["dic"][kind][vid - 1]; nv = NV[it["type"]]
                                    dd = {}; _dom(b, p5, s5, nv, dd)   # sub-chunks (region id, floats); node data uses id 0 -> key -1
                                    st["vars"][it["name"]] = dd[-1] if kind == "node" else dd
            out["states"].append(st)
    out["truncated"] = trunc; return out
def _dom(b, p, s, nv, dct):
    for d, p6, s6 in chunks(b, p, p + s): dct[d - 1] = np.frombuffer(b, np.float32, s6 // 4, p6).reshape(-1, nv).astype(float)

