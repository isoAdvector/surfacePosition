#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Prepares a copy of a case to be run again without function objects: removes
# the functions of system/controlDict, the time directories after 0 (also of
# the processor directories), postProcessing and the solver log.
#
# Usage: withoutFunctionObjects.py <caseDir>

import os
import re
import shutil
import sys

case = sys.argv[1]


def is_later_time(name):
    try:
        return float(name) > 0
    except ValueError:
        return False


for d in [case] + [os.path.join(case, p) for p in os.listdir(case) if p.startswith("processor")]:
    for name in os.listdir(d):
        if is_later_time(name):
            shutil.rmtree(os.path.join(d, name))

shutil.rmtree(os.path.join(case, "postProcessing"), ignore_errors=True)
for name in os.listdir(case):
    if name.startswith("log.") and name not in ("log.blockMesh", "log.setAlphaField", "log.decomposePar"):
        os.remove(os.path.join(case, name))

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
