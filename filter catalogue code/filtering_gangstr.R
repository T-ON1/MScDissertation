library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/GangSTR-master"
bed_file <- file.path(base_dir, "hg38_ver17.bed")

filtered_1_output_file <- file.path(base_dir, "gangstr_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "gangstr_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "gangstr_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "gangstr_removed_2.0.tsv")

catalogue_name <- "gangstr_hg38_ver17"

gangstr_columns <- c(
  "chrom", "start", "end", "period", "motif", "repeat_sequence"
)

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)


gangstr_catalogue <- read_tsv(
  bed_file,
  col_names = gangstr_columns,
  quote = "",
  col_types = cols(
    chrom = col_character(),
    start = col_integer(),
    end = col_integer(),
    period = col_integer(),
    motif = col_character(),
    repeat_sequence = col_character()
  ),
  show_col_types = FALSE,
  progress = TRUE
) %>%
  transmute(
    catalogue = catalogue_name,
    region_chrom = chrom,
    region_start = start,
    region_end = end,
    motif = motif,
    motif_length = period,
    motif_copies = (region_end - region_start + 1) / motif_length,
    repeat_length = region_end - region_start + 1,
    repeat_sequence = repeat_sequence
  )

gangstr_filtered_1 <- gangstr_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

gangstr_removed_1 <- gangstr_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(gangstr_filtered_1, filtered_1_output_file)
write_tsv(gangstr_removed_1, removed_1_output_file)



gangstr_step_2 <- gangstr_filtered_1 %>%
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

gangstr_filtered_2 <- gangstr_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

gangstr_removed_2 <- gangstr_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(gangstr_filtered_2, filtered_2_output_file)
write_tsv(gangstr_removed_2, removed_2_output_file)
