#!/usr/bin/env python3

import pandas as pd
import platform

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


sizes_path = "${sizes}"

df = pd.read_csv(sizes_path, sep="\\t")

df = df.pivot(columns="state", index="sample", values="size")

state_order = ["unfiltered", "filtered", "thresholded", "dedoubleted"]
state_order = [col for col in state_order if col in df.columns]

df = df[state_order].T

# Add a total column
df["total"] = df.sum(axis=1)

df.to_csv("${prefix}.tsv", sep="\\t")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
