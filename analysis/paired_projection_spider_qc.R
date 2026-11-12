# Basic data-admission QC only. No gene-target association, model, or prediction.
# Missing reporter values stay unknown; zero capture is not anatomical absence.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
o <- readRDS(args[[1]])
md <- attr(o, "meta.data")
counts <- attr(attr(o, "assays")$RNA, "counts")
d <- attr(counts, "Dim")
values <- attr(counts, "x")
stopifnot(is.data.frame(md), is.numeric(values), "LHA-I" %in% names(md))
lha <- md[["LHA-I"]]
stopifnot(is.numeric(lha))
qc <- c(
  cell_rows = nrow(md),
  unique_cell_ids = length(unique(rownames(md))),
  gene_rows = d[[1]],
  gene_columns = d[[2]],
  stored_gene_values = length(values),
  nonfinite_stored_gene_values = sum(!is.finite(values)),
  negative_stored_gene_values = sum(values < 0, na.rm = TRUE),
  lha_missing = sum(is.na(lha)),
  lha_nonfinite_nonmissing = sum(!is.finite(lha) & !is.na(lha)),
  lha_negative = sum(lha < 0, na.rm = TRUE),
  lha_positive_capture = sum(lha > 0, na.rm = TRUE),
  lha_zero_capture = sum(lha == 0, na.rm = TRUE),
  slice_labels = length(unique(md$slice)),
  explicit_animal_id_fields = length(grep("animal|mouse|sample|subject|batch", names(md), ignore.case = TRUE))
)
write.table(data.frame(check = names(qc), value = as.numeric(qc)), args[[2]],
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Basic QC saved; no marker-target relationships evaluated.\n")
