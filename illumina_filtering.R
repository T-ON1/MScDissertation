library(tidyverse)
library(jsonlite)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/illumina 174K"
json_file <- file.path(base_dir, "RepeatCatalogs-master", "hg38", "variant_catalog.json")

filtered_1_output_file <- file.path(base_dir, "illumina_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "illumina_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "illumina_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "illumina_removed_2.0.tsv")

catalogue_name <- "illumina_174K_hg38"

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)

extract_motifs <- function(locus_structure) 
  motifs <- str_match_all(locus_structure, "\\(([^)]+)\\)[*+]")[[1]][, 2]

catalogue_raw <- fromJSON(json_file, simplifyVector = FALSE)

illumina_catalogue <- map_dfr(catalogue_raw, function(locus) {
  reference_regions <- as_vector(locus$ReferenceRegion)
  n_regions <- length(reference_regions)
  variant_ids <- get_or_default(locus$VariantId, locus$LocusId)
  variant_types <- get_or_default(locus$VariantType, NA_character_)

  tibble(
    catalogue = catalogue_name,
    locus_id = locus$LocusId,
    variant_id = rep_len(variant_ids, n_regions),
    variant_type = rep_len(variant_types, n_regions),
    locus_structure = locus$LocusStructure,
    reference_region = reference_regions,
    motif = rep_len(extract_motifs(locus$LocusStructure), n_regions)
  )
}) %>%
  extract(
    reference_region,
    into = c("region_chrom", "region_start", "region_end"),
    regex = "^([^:]+):([0-9]+)-([0-9]+)$",
    convert = TRUE,
    remove = FALSE
  ) %>%
  mutate(
    motif_length = nchar(motif),
    repeat_length = region_end - region_start,
    motif_copies = repeat_length / motif_length
  ) %>%
  select(
    catalogue,
    locus_id,
    variant_id,
    variant_type,
    region_chrom,
    region_start,
    region_end,
    motif,
    motif_length,
    motif_copies,
    repeat_length,
    locus_structure,
    reference_region
  )

illumina_filtered_1 <- illumina_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

illumina_removed_1 <- illumina_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(illumina_filtered_1, filtered_1_output_file)
write_tsv(illumina_removed_1, removed_1_output_file)


illumina_step_2 <- illumina_filtered_1 %>%
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

illumina_filtered_2 <- illumina_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

illumina_removed_2 <- illumina_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(illumina_filtered_2, filtered_2_output_file)
write_tsv(illumina_removed_2, removed_2_output_file)

