#!/usr/bin/env python3
#
# Writes the README.md of the repository, a Markdown version of the report,
# from the template README.md.in and generated/results.json (makeData.py):
#   {{name}}         a number or text of the report (a macro of results.tex)
#   {{rows:name}}    the rows of a table of the report
#
# Usage: makeReadme.py   (from the report directory, after makeData.py)

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
README = os.path.abspath(os.path.join(HERE, "..", "..", "README.md"))

data = json.load(open(os.path.join(HERE, "generated", "results.json")))
text = open(os.path.join(HERE, "README.md.in")).read()

missing = []


def fill(m):
    key = m.group(1)
    if key.startswith("rows:"):
        rows = data["tables"].get(key[5:])
        if rows is None:
            missing.append(key)
            return m.group(0)
        return "\n".join(rows)
    if key not in data["macros"]:
        missing.append(key)
        return m.group(0)
    return data["macros"][key]


text = re.sub(r"\{\{([\w:]+)\}\}", fill, text)
if missing:
    sys.exit(f"Unknown placeholders in README.md.in: {', '.join(missing)}")

with open(README, "w") as f:
    f.write(text)
print(f"Wrote {README}")
