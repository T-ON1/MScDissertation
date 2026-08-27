library(tidyverse)

# This is for the UCSC simpleRepeat table.
# This file has no header row, so the UCSC simpleRepeat column names are added here.
# UCSC simpleRepeat does not contain gene/exon/feature annotation columns.
# The raw name column is retained as ucsc_name.
# UCSC period is retained as source_motif_length, but motif_length is calculated
# from the actual motif sequence so filtering removes motifs longer than 6 bp.
# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/UCSC"
simple_repeat_file <- file.path(base_dir, "simpleRepeat.txt")

filtered_1_output_file <- file.path(base_dir, "ucsc_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "ucsc_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "ucsc_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "ucsc_removed_2.0.tsv")

catalogue_name <- "ucsc_simpleRepeat"


simple_repeat_columns <- c(
  "bin", "chrom", "chromStart", "chromEnd", "name", "period", "copyNum",
  "consensusSize", "perMatch", "perIndel", "score", "A", "C", "G", "T",
  "entropy", "sequence"
)

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)


ucsc_catalogue <- read_tsv(
  simple_repeat_file,
  col_names = simple_repeat_columns,
  quote = "",
  col_types = cols(
    chrom = col_character(),
    chromStart = col_integer(),
    chromEnd = col_integer(),
    name = col_character(),
    period = col_double(),
    copyNum = col_double(),
    sequence = col_character(),
    .default = col_skip()
   )
) %>%
  transmute(
    catalogue = catalogue_name,
    ucsc_name = name,
    region_chrom = chrom,
    region_start = chromStart,
    region_end = chromEnd,
    motif = sequence,
    source_motif_length = as.integer(period),
    motif_length = nchar(motif),
    motif_copies = copyNum,
    calculated_motif_copies = (region_end - region_start) / motif_length,
    repeat_length = region_end - region_start
  )

ucsc_filtered_1 <- ucsc_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

ucsc_removed_1 <- ucsc_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(ucsc_filtered_1, filtered_1_output_file)
  write_tsv(ucsc_removed_1, removed_1_output_file)


ucsc_step_2 <- ucsc_filtered_1 %>%
  mutate(
    minimum_repeat_length = case_when(
      motif_length == 1L ~ 9L,
      motif_length == 2L ~ 10L,
      motif_length == 3L ~ 9L,
      motif_length == 4L ~ 12L,
      motif_length == 5L ~ 15L,
      motif_length == 6L ~ 18L,
      TRUE ~ NA_integer_
    )
  )

ucsc_filtered_2 <- ucsc_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

ucsc_removed_2 <- ucsc_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(ucsc_filtered_2, filtered_2_output_file)
write_tsv(ucsc_removed_2, removed_2_output_file)