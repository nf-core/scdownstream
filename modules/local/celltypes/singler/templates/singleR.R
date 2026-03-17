#!/usr/bin/env Rscript

library(SingleR)
library(celldex)
library(yaml)
library(Seurat)
library(ggplot2)
library(anndataR)
library(HDF5Array)

# Function to convert a list to a YAML-like structure, as before
format_yaml_like <- function(data, indent = 0) {
  yaml_str <- ""
  for (key in names(data)) {
    spaces <- strrep("  ", indent)
    value <- data[[key]]
    if (is.list(value)) {
      yaml_str <- paste0(yaml_str, spaces, key, ":\n",
                        format_yaml_like(value, indent + 1))
    } else {
      yaml_str <- paste0(yaml_str, spaces, key, ": ", value, "\n")
    }
  }
  return(yaml_str)
}
# Read .h5ad file using zellkonverter
symbol_col <- "${symbol_col}"
h5ad_file <- "${h5ad}" # Get the filename from environment variable
sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment") # Converts .h5ad to a SingleCellExperiment object

if (symbol_col != "index") {
  if (!symbol_col %in% colnames(rowData(sce))) {
    stop(paste0("Symbol column ", symbol_col, " not found in adata.var.columns"))
  }
  rownames(sce) <- rowData(sce)[[symbol_col]]
}

# Split the references by comma and loop over each
prefix <- "${prefix}"
num_threads <- max(1L, as.integer("${task.cpus}"))
references <- strsplit("${references.join(',')}", ",")[[1]]
reference_names <- strsplit("${names.join(',')}", ",")[[1]]
reference_labels <- strsplit("${labels.join(',')}", ",")[[1]]

stopifnot(
    #"Lengths of references and reference_labels vectors must match",
    length(references) == length(reference_labels) && length(references) == length(reference_names)
)

Sys.setenv(XDG_CACHE_HOME = file.path(getwd(), ".cache"))
prediction_results <- list()
for (ref_idx in seq_along(references)) {
  ref <- references[ref_idx]
  reflabel <- reference_labels[ref_idx]
  ref_name <- reference_names[ref_idx]

  # Untar the reference files into a directory named after the reference without the extension
  ref_dir <- gsub(".tar", "", ref)
  untar(ref, exdir = "./")
  # Read the SummarizedExperiment object from the provided path
  reference <- loadHDF5SummarizedExperiment(dir = ref_dir)
  stopifnot(
    reflabel %in% colnames(colData(reference))
  )
  predictions <- SingleR(
    test = assay(sce, 'counts'),
    ref = reference,
    labels = colData(reference)[[reflabel]],
    num.threads = num_threads
  )

  # Plot and save heatmap
  p <- plotScoreHeatmap(
    predictions,
    main = paste0(
      "SingleR Predictions: ",
      basename(h5ad_file),
      " [", prefix, "_", ref_name, "]"
    ),
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  ggsave(
    filename = paste0(prefix, "_", ref_name, "_heatmap.pdf"),
    plot = p,
    width = 10,
    height = 8
  )

  # Plot and save distribution
  p2 <- plotDeltaDistribution(predictions, ncol = 3)
  p2 <- p2 + ggtitle(
    paste0(
      "SingleR Predictions: ",
      basename(h5ad_file),
      " [", prefix, "_", ref_name, "]"
    )
  )
  ggsave(
    filename = paste0(prefix, "_", ref_name, "_distribution.pdf"),
    plot = p2,
    width = 14,
    height = 12
  )

  colnames(predictions) <- paste0(
    colnames(predictions), "_", ref_name
  )
  prediction_results[[ref]] <- predictions
}

prediction_nrows <- lapply(prediction_results, nrow)
prediction_rownames <- lapply(prediction_results, rownames)


stopifnot(
  all(sapply(prediction_nrows, function(x) x == prediction_nrows[[1]])) &
  all(sapply(prediction_rownames, function(x) all(x == prediction_rownames[[1]])))
)

# This is predicated in the assumption that all prediction data frames have exactly
# the same rows ... see the stopifnot clause above
predictions <- do.call(cbind, prediction_results)

write.csv(
  predictions,
  file = paste0(prefix, "_predictions.csv"),
  row.names = TRUE
)

# Capturing version information, as before
versions <- list(
  "${task.process}" = list(
    R = R.version.string,
    SingleR = as.character(packageVersion("SingleR")),
    celldex = as.character(packageVersion("celldex")),
    anndataR = as.character(packageVersion("anndataR")),
    yaml = as.character(packageVersion("yaml")),
    Seurat = as.character(packageVersion("Seurat")),
    ggplot2 = as.character(packageVersion("ggplot2"))
  )
)

# Delete the Rplots.pdf file if it exists
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# Write versions info into a YAML file, as before
write(format_yaml_like(versions), file = "versions.yml")
