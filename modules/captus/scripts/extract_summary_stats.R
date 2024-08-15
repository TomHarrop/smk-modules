#!/usr/bin/env Rscript


log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


get_summary_stats <- function(
    extraction_stats,
    id_cols = id_cols,
    summary_cols = summary_cols) {
  sample_medians <- extraction_stats[
    , lapply(.SD, median, na.rm = TRUE),
    by = id_cols,
    .SDcols = summary_cols
  ]
  sample_means <- extraction_stats[
    , lapply(.SD, mean, na.rm = TRUE),
    by = id_cols,
    .SDcols = summary_cols
  ]

  summary_stats <- merge(
    sample_medians,
    sample_means,
    by = id_cols,
    suffixes = c("_median", "_mean")
  )
  setkey(
    summary_stats,
    marker_type,
    sample_name
  )
  return(summary_stats)
}

###########
# GLOBALS #
###########


# INPUTS #
stats_file <- snakemake@params[["captus_extraction_stats"]]
# OUTPUTS #
sample_stats_by_pct_recovered_file <-
  snakemake@output[["sample_stats_by_pct_recovered"]]
sample_stats_by_wscore_file <- snakemake@output[["sample_stats_by_wscore"]]

# summary stats
id_cols <- c(
  "sample_name",
  "marker_type"
)

summary_cols <- c(
  "ref_len_matched",
  "pct_recovered",
  "pct_identity",
  "score",
  "wscore",
  "hit_len",
  "cds_len",
  "intron_len",
  "flanks_len"
)

marker_order <- c(
  "NUC", "PTD", "MIT", "CLR"
)

########
# MAIN #
########

library(data.table)

captus_extraction_stats <- fread(stats_file)

# order the output
captus_extraction_stats[
  ,
  marker_type := factor(marker_type, levels = marker_order)
]

# only count one reference per locus
# ad-hoc method using pct_recovered
best_hits_by_pct_recovered <- captus_extraction_stats[
  , .SD[which.max(pct_recovered)],
  by = c(id_cols, "locus")
]
sample_stats_by_pct_recovered <-
  get_summary_stats(best_hits_by_pct_recovered, id_cols, summary_cols)


# Captus's method using pct_recovered and identity,
# https://edgardomortiz.github.io/captus.docs/assembly/extract/output/#fasta-headers-explanation
best_hits_by_captus_wscore <- captus_extraction_stats[
  , .SD[hit == 0],
  by = c(id_cols, "locus")
]
sample_stats_by_wscore <-
  get_summary_stats(best_hits_by_captus_wscore, id_cols, summary_cols)

# write the output
fwrite(sample_stats_by_pct_recovered, sample_stats_by_pct_recovered_file)
fwrite(sample_stats_by_wscore, sample_stats_by_wscore_file)
sessionInfo()
