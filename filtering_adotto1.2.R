library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/adottoTR_filtering"
bed_file <- file.path(base_dir, "adotto_TRregions_v1.2.bed")

filtered_1_output_file <- file.path(base_dir, "adotto_TRregions_v1.2_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "adotto_TRregions_v1.2_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "adotto_TRregions_v1.2_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "adotto_TRregions_v1.2_removed_2.0.tsv")

catalogue_name <- "adotto_TRregions_v1.2"


extract_matches <- function(x, pattern, numeric = FALSE) {
  matches <- str_match_all(x, pattern)

  map(matches, function(match) {
    values <- match[, 2]

    if (numeric) {
      as.numeric(values)
    } else {
      values
    }
  })
}

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)

adotto_raw <- read_tsv(
  bed_file,
  col_names = paste0("X", 1:18),
  quote = "\"",
  col_types = cols(
    X13 = col_character(),
    X17 = col_character(),
    X18 = col_character(),
    .default = col_skip()
  )
)

adotto_catalogue <- adotto_raw %>%
  transmute(
    catalogue = catalogue_name,
    repeat_class = na_if(X13, "."),
    genomic_feature = na_if(X17, "."),
    region_chrom = extract_matches(X18, '"chrom"\\s*:\\s*"([^"]+)"'),
    region_start = extract_matches(X18, '"start"\\s*:\\s*([0-9]+)', numeric = TRUE),
    region_end = extract_matches(X18, '"end"\\s*:\\s*([0-9]+)', numeric = TRUE),
    motif = extract_matches(X18, '"motif"\\s*:\\s*"([^"]+)"'),
    motif_copies = extract_matches(X18, '"copies"\\s*:\\s*([0-9.]+)', numeric = TRUE)
    ) %>%
  unnest(c(region_chrom, region_start, region_end, motif, motif_copies)) %>%
  mutate(
    motif_length = nchar(motif),
    repeat_length = (region_end - region_start) + 1,
    calculated_motif_copies = repeat_length / motif_length
  ) %>%
  select(
    catalogue,
    repeat_class,
    genomic_feature,
    region_chrom,
    region_start,
    region_end,
    motif,
    motif_length,
    motif_copies,
    calculated_motif_copies,
    repeat_length
)

adotto_filtered_1 <- adotto_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

adotto_removed_1 <- adotto_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(adotto_filtered_1, filtered_1_output_file)
write_tsv(adotto_removed_1, removed_1_output_file)



adotto_step_2 <- adotto_filtered_1 %>%
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

adotto_filtered_2 <- adotto_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

adotto_removed_2 <- adotto_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(adotto_filtered_2, filtered_2_output_file)
write_tsv(adotto_removed_2, removed_2_output_file)
