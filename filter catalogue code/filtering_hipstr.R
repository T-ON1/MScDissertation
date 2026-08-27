library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/HipSTR"
bed_file <- file.path(base_dir, "hg38.hipstr_reference.bed")

filtered_1_output_file <- file.path(base_dir, "hipstr_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "hipstr_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "hipstr_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "hipstr_removed_2.0.tsv")

catalogue_name <- "hipstr_hg38"

hipstr_columns <- c(
  "chrom", "start", "end", "period", "copies", "locus_id", "motif"
)

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)
required_hipstr_output_columns <- c("locus_id", "source_motif_length")


hipstr_catalogue <- read_tsv(
  bed_file,
  col_names = hipstr_columns,
  quote = "",
  col_types = cols(
    chrom = col_character(),
    start = col_integer(),
    end = col_integer(),
    period = col_integer(),
    copies = col_double(),
    locus_id = col_character(),
    motif = col_character()
  ),
  show_col_types = FALSE,
  progress = TRUE
) %>%
  mutate(
    motif = na_if(motif, "N/A"),
    motif_parts = str_split(motif, fixed("/")),
    motif_length = map_int(
      motif_parts,
      function(parts) {
        if (any(is.na(parts))) {
          return(NA_integer_)
        }

        part_lengths <- nchar(parts)

        if (n_distinct(part_lengths) == 1L) {
          return(part_lengths[[1]])
        }

        NA_integer_
      }
    )
  ) %>%
  transmute(
    catalogue = catalogue_name,
    locus_id = locus_id,
    region_chrom = chrom,
    region_start = start,
    region_end = end,
    motif = motif,
    motif_length = motif_length,
    source_motif_length = period,
    motif_copies = copies,
    repeat_length = region_end - region_start + 1,
    calculated_motif_copies = repeat_length / motif_length
)

hipstr_filtered_1 <- hipstr_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

hipstr_removed_1 <- hipstr_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(hipstr_filtered_1, filtered_1_output_file)
write_tsv(hipstr_removed_1, removed_1_output_file)

hipstr_step_2 <- hipstr_filtered_1 %>%
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

hipstr_filtered_2 <- hipstr_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

hipstr_removed_2 <- hipstr_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(hipstr_filtered_2, filtered_2_output_file)
write_tsv(hipstr_removed_2, removed_2_output_file)
