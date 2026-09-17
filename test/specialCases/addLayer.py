#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Adds a layer of water above the water of the layers scenario, in
# 0/alpha.water of the 10x10x10 unit cube (cell label i + 10 j + 100 k):
# cells with 0.6 < z < 0.7 full and 0.7 < z < 0.8 half full, so that the
# layer spans 0.6 < z < 0.75 with a sharp lower interface.
#
# Usage: addLayer.py   (in the case directory, after setAlphaField)

import re

path = "0/alpha.water"
text = open(path).read()
m = re.search(r"(internalField\s+nonuniform\s+List<scalar>\s*)(\d+)(\s*\()([^)]*)(\))", text)
values = m.group(4).split()
assert len(values) == 1000, len(values)

for k, value in ((6, "1"), (7, "0.5")):
    for j in range(10):
        for i in range(10):
            values[i + 10*j + 100*k] = value

text = text[:m.start(4)] + "\n" + "\n".join(values) + "\n" + text[m.end(4):]
open(path, "w").write(text)
