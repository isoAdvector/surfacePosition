#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Removes the function objects of a case: empties the functions of
# system/controlDict.
#
# Usage: withoutFunctionObjects.py <caseDir>

import os
import re
import sys

case = sys.argv[1]

path = os.path.join(case, "system", "controlDict")
text = open(path).read()
m = re.search(r"^functions\s*\{", text, re.M)
if m:
    depth = 0
    for i in range(m.end() - 1, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                text = text[:m.start()] + "functions\n{\n}" + text[i + 1:]
                break
open(path, "w").write(text)
