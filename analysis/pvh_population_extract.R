# Own base-R reader for the RNA "counts" layer of one author Seurat object (PB or SC).
# No author code is sourced; no Seurat/SeuratObject/Matrix functions are called; no package is installed.
# Structure is read with attr() only. Only the requested cohort cells are exported.
#
# Usage: Rscript --vanilla analysis/pvh_population_extract.R <rds> <object PB|SC> <expected_md5> <cohort.tsv> <out_dir>
#   cohort.tsv: two tab-separated columns with header "object<TAB>cell_id"
#   out_dir must not exist. Writes genes.txt, cells.tsv, triplets.tsv, receipt.tsv.

refuse <- function(...) {
  message("REFUSED: ", ...)
  quit(save = "no", status = 3L)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) refuse("expected 5 arguments")
rds <- args[[1L]]; object_name <- args[[2L]]; expected_md5 <- args[[3L]]
cohort_path <- args[[4L]]; out_dir <- args[[5L]]
if (!(object_name %in% c("PB", "SC"))) refuse("object must be PB or SC")
if (!file.exists(rds)) refuse("rds not found")
if (file.exists(out_dir)) refuse("output directory already exists")
if (unname(tools::md5sum(rds)) != expected_md5) refuse("rds md5 differs from the pin")

cohort <- read.delim(cohort_path, colClasses = "character", quote = "", comment.char = "")
if (!identical(colnames(cohort), c("object", "cell_id"))) refuse("cohort columns must be object, cell_id")
wanted <- cohort$cell_id[cohort$object == object_name]
if (length(wanted) < 1L) refuse("no cohort cells for this object")
if (anyDuplicated(wanted)) refuse("duplicate cohort cell IDs")

object <- readRDS(rds)
meta <- attr(object, "meta.data", exact = TRUE)
if (!is.data.frame(meta)) refuse("meta.data missing")
meta_cells <- rownames(meta)
if (is.null(meta_cells) || anyDuplicated(meta_cells)) refuse("meta.data row names missing or duplicated")

assays <- attr(object, "assays", exact = TRUE)
if (!("RNA" %in% names(assays))) refuse("RNA assay missing")
assay <- assays[["RNA"]]
layers <- attr(assay, "layers", exact = TRUE)
if (!("counts" %in% names(layers))) refuse("counts layer missing")
layer <- layers[["counts"]]

i <- attr(layer, "i", exact = TRUE); p <- attr(layer, "p", exact = TRUE)
x <- attr(layer, "x", exact = TRUE); dims <- attr(layer, "Dim", exact = TRUE)
if (is.null(i) || is.null(p) || is.null(x) || length(dims) != 2L) refuse("counts layer is not a column-compressed sparse matrix")
n_genes <- dims[[1L]]; n_cells <- dims[[2L]]
if (length(p) != n_cells + 1L || p[[length(p)]] != length(x) || length(i) != length(x)) refuse("inconsistent sparse structure")
if (any(!is.finite(x)) || any(x < 0) || any(x != round(x))) refuse("counts layer is not non-negative whole numbers")
if ("data" %in% names(layers)) {
  dx <- attr(layers[["data"]], "x", exact = TRUE)
  if (!is.null(dx) && identical(as.numeric(dx), as.numeric(x))) refuse("counts layer identical to data layer")
}

map_names <- function(map, layer_name) {
  if (is.null(map)) return(NULL)
  m <- tryCatch(unclass(map), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  dn <- attr(m, "dimnames", exact = TRUE); dm <- attr(m, "dim", exact = TRUE)
  if (is.null(dn) || length(dm) != 2L || !(layer_name %in% dn[[2L]])) return(NULL)
  v <- as.logical(as.vector(m))
  if (length(v) != dm[[1L]] * dm[[2L]] || any(is.na(v))) return(NULL)
  col <- v[(match(layer_name, dn[[2L]]) - 1L) * dm[[1L]] + seq_len(dm[[1L]])]
  dn[[1L]][col]
}
layer_dn <- attr(layer, "Dimnames", exact = TRUE)
genes <- map_names(attr(assay, "features", exact = TRUE), "counts")
if (is.null(genes) && !is.null(layer_dn)) genes <- layer_dn[[1L]]
cells <- map_names(attr(assay, "cells", exact = TRUE), "counts")
if (is.null(cells) && !is.null(layer_dn)) cells <- layer_dn[[2L]]
if (is.null(genes) || length(genes) != n_genes || anyDuplicated(genes)) refuse("gene IDs for the counts layer could not be established exactly")
if (is.null(cells) || length(cells) != n_cells || !identical(as.character(cells), as.character(meta_cells))) {
  refuse("counts-layer cell order could not be matched exactly to meta.data row names")
}

columns <- match(wanted, cells)
if (any(is.na(columns))) refuse("cohort cell not found in the counts layer: ", wanted[is.na(columns)][[1L]])

dir.create(out_dir, recursive = FALSE)
writeLines(genes, file.path(out_dir, "genes.txt"))
write.table(data.frame(cell_index = seq_along(wanted), cell_id = wanted, layer_column = columns),
            file.path(out_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lens <- p[columns + 1L] - p[columns]
idx <- unlist(lapply(columns, function(j) if (p[[j + 1L]] > p[[j]]) (p[[j]] + 1L):p[[j + 1L]] else integer(0)))
gi <- i[idx] + 1L; ci <- rep.int(seq_along(columns), lens); cv <- x[idx]
write.table(data.frame(gene_index = gi, cell_index = ci, count = format(cv, scientific = FALSE, trim = TRUE)),
            file.path(out_dir, "triplets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
receipt <- c(
  paste("object", object_name, sep = "\t"),
  paste("rds_md5", expected_md5, sep = "\t"),
  paste("assay", "RNA", sep = "\t"),
  paste("layer", "counts", sep = "\t"),
  paste("layer_genes", n_genes, sep = "\t"),
  paste("layer_cells", n_cells, sep = "\t"),
  paste("layer_nonzero", length(x), sep = "\t"),
  paste("exported_cells", length(columns), sep = "\t"),
  paste("exported_nonzero", length(cv), sep = "\t"),
  paste("counts_whole_nonnegative", "TRUE", sep = "\t"),
  paste("R_version", R.version.string, sep = "\t"),
  paste("packages_loaded_by_script", "none", sep = "\t")
)
writeLines(receipt, file.path(out_dir, "receipt.tsv"))
message("extracted ", length(columns), " cells x ", n_genes, " genes from ", object_name)
