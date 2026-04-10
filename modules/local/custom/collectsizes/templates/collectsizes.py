#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import json
import platform

import pandas as pd
import yaml


sizes_path = "${sizes}"

df = pd.read_csv(sizes_path, sep="\\t")

df = df.pivot(columns="state", index="sample", values="size")

state_order = ["unfiltered", "filtered", "thresholded", "dedoubleted"]
state_order = [col for col in state_order if col in df.columns]

df = df[state_order].T

# Add a total column
df["total"] = df.sum(axis=1)

df.to_csv("${prefix}.tsv", sep="\\t")

# MultiQC

with open("${prefix}_mqc.json", "w") as f_json:
    json.dump({
        "id": "${prefix}",
        "plot_type": "table",
        "section_name": "Number of cells",
        "description": "The number of cells present in each sample at different pipeline stages.",
        "data": df.to_dict()
    },
    f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
