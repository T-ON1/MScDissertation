library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.
# The long left/right flank sequence columns are skipped to keep outputs compact.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/popSTR 5,401,401"
popstr_file <- file.path(base_dir, "popSTR_markerInfo_chr1-22_combined.txt")

filtered_1_output_file <- file.path(base_dir, "popstr_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "popstr_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "popstr_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "popstr_removed_2.0.tsv")

catalogue_name <- "popstr_chr1_22"


popstr_columns <- c(
  "chrom",
  "start",
  "end",
  "motif",
  "ref_repeat_count",
  "left_flank_1000bp",
  "right_flank_1000bp",
  "repeat_sequence",
  "min_flank_left",
  "min_flank_right",
  "repeat_purity",
  "marker_slippage",
  "marker_stutter",
  "motif_fraction_A",
  "motif_fraction_C",
  "motif_fraction_G",
  "motif_fraction_T"
)

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)


popstr_catalogue <- read_table(
  popstr_file,
  col_names = popstr_columns,
  col_types = cols(
    chrom = col_character(),
    start = col_integer(),
    end = col_integer(),
    motif = col_character(),
    ref_repeat_count = col_double(),
    left_flank_1000bp = col_skip(),
    right_flank_1000bp = col_skip(),
    repeat_sequence = col_character(),
    min_flank_left = col_integer(),
    min_flank_right = col_integer(),
    repeat_purity = col_double(),
    marker_slippage = col_double(),
    marker_stutter = col_double(),
    motif_fraction_A = col_double(),
    motif_fraction_C = col_double(),
    motif_fraction_G = col_double(),
    motif_fraction_T = col_double()
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
    motif_length = nchar(motif),
    motif_copies = ref_repeat_count,
    calculated_motif_copies = (region_end - region_start + 1) / motif_length,
    repeat_length = region_end - region_start + 1,
    repeat_sequence = repeat_sequence,
    min_flank_left = min_flank_left,
    min_flank_right = min_flank_right,
    repeat_purity = repeat_purity,
    marker_slippage = marker_slippage,
    marker_stutter = marker_stutter,
    motif_fraction_A = motif_fraction_A,
    motif_fraction_C = motif_fraction_C,
    motif_fraction_G = motif_fraction_G,
    motif_fraction_T = motif_fraction_T
  )

popstr_filtered_1 <- popstr_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

popstr_removed_1 <- popstr_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(popstr_filtered_1, filtered_1_output_file)
write_tsv(popstr_removed_1, removed_1_output_file)


popstr_step_2 <- popstr_filtered_1 %>%
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

popstr_filtered_2 <- popstr_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

popstr_removed_2 <- popstr_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(popstr_filtered_2, filtered_2_output_file)
write_tsv(popstr_removed_2, removed_2_output_file)

