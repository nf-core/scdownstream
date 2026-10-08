#!/usr/bin/env python3

import platform
from pathlib import Path

import decoupler as dc
import pandas as pd
import yaml

resource = "${meta.id}"
custom_network = "${custom_network}"
species = "${species}"
prefix = "${prefix}"

OMNIPATH_RESOURCES = ("progeny", "collectri", "hallmark")
# Top PROGENy targets per pathway recommended by the decoupler tutorials
PROGENY_TOP = 500


def normalise_species(value: str) -> str:
    key = value.strip().lower()
    if key in ("homo_sapiens", "hsapiens", "hs", "human", ""):
        return "human"
    if key in ("mus_musculus", "mmusculus", "mm", "mouse"):
        return "mouse"
    return key


def read_gmt(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 3:
                continue
            rows.extend((fields[0], gene) for gene in fields[2:] if gene)
    return pd.DataFrame(rows, columns=["source", "target"])


def read_table(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\\t"
    table = pd.read_csv(path, sep=sep)
    missing = {"source", "target"} - set(table.columns)
    if missing:
        raise ValueError(f"Custom network {path.name} lacks columns: {sorted(missing)}")
    return table


def load_omnipath(name: str, organism: str) -> pd.DataFrame:
    if name == "progeny":
        return dc.op.progeny(organism=organism, top=PROGENY_TOP)
    if name == "collectri":
        return dc.op.collectri(organism=organism)
    if name == "hallmark":
        return dc.op.hallmark(organism=organism)
    raise ValueError(f"Unknown decoupler resource '{name}'; expected one of {OMNIPATH_RESOURCES} or a custom network")


if custom_network:
    path = Path(custom_network)
    net = read_gmt(path) if path.suffix.lower() == ".gmt" else read_table(path)
else:
    net = load_omnipath(resource, normalise_species(species))

if "weight" not in net.columns:
    net["weight"] = 1.0

net = net[["source", "target", "weight"]].dropna(subset=["source", "target"])
net["source"] = net["source"].astype(str)
net["target"] = net["target"].astype(str)
net["weight"] = pd.to_numeric(net["weight"], errors="coerce").fillna(1.0)
net = net.drop_duplicates(["source", "target"]).sort_values(["source", "target"]).reset_index(drop=True)

if net.empty:
    raise ValueError(f"Network '{resource}' contains no interactions")

net.to_csv(f"{prefix}.tsv", sep="\\t", index=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "decoupler": dc.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
