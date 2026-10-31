#!/usr/bin/env python3
"""subset_json.py IN.json OUT.json key [key ...]: writes the listed top-level keys."""
import sys, json
d = json.load(open(sys.argv[1])); json.dump({k: d[k] for k in sys.argv[3:]}, open(sys.argv[2], "w"))
