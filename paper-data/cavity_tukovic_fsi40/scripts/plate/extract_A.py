"""Print pointD at the top mid-span vertex (x=4, y=0, z=0) of a plate run.
Usage: python3 -I extract_A.py CASEDIR NX NY [TIME]"""
import sys, re
case, nx, ny = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
t = sys.argv[4] if len(sys.argv) > 4 else "1"
txt = open(f"{case}/{t}/pointD").read()
m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)\s*\n?\(", txt)
n = int(m.group(1)); body = txt[m.end():]
vals = re.findall(r"\(([^()]*)\)", body)[:n]
i = nx // 2 + ny * (nx + 1)
print(case, vals[i])
