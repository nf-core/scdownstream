#!/usr/bin/env python3

import importlib.metadata
import platform

import anndata as ad
import numpy as np
import pandas as pd
import yaml

LNCRNA_BIOTYPES = {
    "lncRNA",
    "lincRNA",
    "antisense",
    "antisense_RNA",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "3prime_overlapping_ncRNA",
    "3prime_overlapping_ncrna",
    "bidirectional_promoter_lncRNA",
    "macro_lncRNA",
    "non_coding",
}
SMALL_NCRNA_BIOTYPES = {
    "miRNA",
    "snRNA",
    "snoRNA",
    "scaRNA",
    "misc_RNA",
    "sRNA",
    "scRNA",
    "vault_RNA",
    "vaultRNA",
    "Y_RNA",
    "ribozyme",
    "rRNA",
    "tRNA",
}
COLUMNS = ["gene_biotype", "chromosome"]


def strip_version(gene_id):
    head, sep, tail = gene_id.rpartition(".")
    return head if sep and tail.isdigit() else gene_id


def normalise_chromosome(value):
    value = str(value)
    if value.lower().startswith("chr"):
        value = value[3:]
    return "MT" if value.upper() in {"M", "MT"} else value


def is_secondary_contig(chromosome):
    # Patches, alternate haplotypes and unplaced scaffolds (e.g. HSCHR6_MHC_COX_CTG1, KI270728.1)
    return isinstance(chromosome, str) and ("_" in chromosome or "." in chromosome)


def biotype_group(biotype, chromosome):
    if chromosome == "MT" or (isinstance(biotype, str) and biotype.lower().startswith("mt_")):
        return "mitochondrial"
    if not isinstance(biotype, str):
        return np.nan
    if "pseudogene" in biotype:
        return "pseudogene"
    if biotype == "protein_coding":
        return "protein_coding"
    if biotype.startswith(("IG_", "TR_")) and biotype.endswith("_gene"):
        return "IG_TR_gene"
    if biotype in LNCRNA_BIOTYPES:
        return "lncRNA"
    if biotype in SMALL_NCRNA_BIOTYPES:
        return "small_ncRNA"
    return "other"


def candidates_from_table(path, genes):
    table = pd.read_csv(path, sep="\\t", dtype=str)
    biotype_col = next((col for col in ["gene_biotype", "gene_type"] if col in table.columns), None)
    missing = [col for col in ["gene_id", "chromosome"] if col not in table.columns]
    if biotype_col is None or missing:
        raise ValueError(
            "Gene table must contain gene_id, chromosome and gene_biotype or gene_type columns. "
            f"Found: {', '.join(table.columns)}"
        )

    records = pd.DataFrame(
        {
            "gene_id": table["gene_id"].map(strip_version),
            "gene_name": table["gene_name"] if "gene_name" in table.columns else np.nan,
            "gene_biotype": table[biotype_col],
            "chromosome": table["chromosome"],
        }
    )
    queries = pd.DataFrame({"query": genes, "gene_id": [strip_version(gene) for gene in genes]})

    # Match on gene name first; genes without a name match fall back to the gene ID.
    by_name = queries[["query"]].merge(records.dropna(subset=["gene_name"]), left_on="query", right_on="gene_name")
    unmatched = queries[~queries["query"].isin(by_name["query"])]
    by_id = unmatched.merge(records, on="gene_id")
    return pd.concat([by_name, by_id], ignore_index=True)[["query"] + COLUMNS]


def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def candidates_from_mygene(genes, species):
    import mygene

    hits = mygene.MyGeneInfo().querymany(
        list(genes),
        scopes=["symbol", "ensembl.gene"],
        fields="symbol,ensembl.type_of_gene,genomic_pos.chr",
        species=species,
        verbose=False,
    )
    hits = [hit for hit in hits if not hit.get("notfound")]

    # A symbol query also returns genes that only list it as an alias or previous symbol.
    # Keep exact symbol matches when there are any; ID queries keep all of their hits.
    exact = {hit["query"] for hit in hits if str(hit.get("symbol", "")).casefold() == str(hit["query"]).casefold()}
    hits = [
        hit
        for hit in hits
        if hit["query"] not in exact or str(hit.get("symbol", "")).casefold() == hit["query"].casefold()
    ]

    rows = []
    for hit in hits:
        biotypes = [entry.get("type_of_gene") for entry in as_list(hit.get("ensembl")) if isinstance(entry, dict)]
        chromosomes = [entry.get("chr") for entry in as_list(hit.get("genomic_pos")) if isinstance(entry, dict)]
        for biotype in biotypes or [None]:
            for chromosome in chromosomes or [None]:
                rows.append({"query": hit["query"], "gene_biotype": biotype, "chromosome": chromosome})
    return pd.DataFrame(rows, columns=["query"] + COLUMNS)


def summarise(candidates, genes):
    candidates = candidates.dropna(subset=COLUMNS, how="all").copy()
    candidates["chromosome"] = candidates["chromosome"].map(normalise_chromosome, na_action="ignore")

    # Ignore secondary contigs for genes that also have a primary chromosome entry
    secondary = candidates["chromosome"].map(is_secondary_contig).astype(bool)
    has_primary = candidates["chromosome"].notna() & ~secondary
    gene_has_primary = has_primary.groupby(candidates["query"]).transform("any")
    candidates = candidates[~(secondary & gene_has_primary)]

    result = pd.DataFrame(index=pd.Index(genes))
    ambiguous = pd.Series(False, index=result.index)
    for column in COLUMNS:
        values = candidates.dropna(subset=[column]).drop_duplicates(["query", column])
        counts = values.groupby("query")[column].size()
        unique = values[values["query"].map(counts) == 1].set_index("query")[column]
        result[column] = unique.reindex(result.index)
        ambiguous |= (counts > 1).reindex(result.index, fill_value=False)

    result["biotype_group"] = [
        biotype_group(biotype, chromosome)
        for biotype, chromosome in zip(result["gene_biotype"], result["chromosome"], strict=True)
    ]
    result["gene_annotation_ambiguous"] = ambiguous.astype(bool)
    return result


adata = ad.read_h5ad("${h5ad}", backed="r")
genes = adata.var_names.astype(str).tolist()
if len(set(genes)) != len(genes):
    raise ValueError("The var index contains duplicate gene names")

gene_table = "${gene_table}"
if gene_table:
    candidates = candidates_from_table(gene_table, genes)
else:
    candidates = candidates_from_mygene(genes, "${species}")

result = summarise(candidates, genes)
result.to_csv("${prefix}.csv")

# The pickle keeps dtypes for ADATA_EXTEND; a CSV would let pandas re-infer e.g. chromosome as numeric
for column in ["gene_biotype", "chromosome", "biotype_group"]:
    result[column] = result[column].astype("category")
result.to_pickle("${prefix}.pkl")

n_annotated = int(result["gene_biotype"].notna().sum())
print(f"Annotated {n_annotated} of {len(genes)} genes; {int(result['gene_annotation_ambiguous'].sum())} ambiguous")
print(result["biotype_group"].value_counts(dropna=False).to_string())

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": importlib.metadata.version("pandas"),
    }
}
if not gene_table:
    versions["${task.process}"]["mygene"] = importlib.metadata.version("mygene")

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
