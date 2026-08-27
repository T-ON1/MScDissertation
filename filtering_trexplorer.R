library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat/tract lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/TRexplorer"
bed_file <- file.path(base_dir, "trexplorer_bed_layout.bed")

filtered_1_output_file <- file.path(base_dir, "trexplorer_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "trexplorer_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "trexplorer_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "trexplorer_removed_2.0.tsv")

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)

required_annotation_columns <- c(
  "locus_id",
  "variant_type",
  "source",
  "canonical_motif",
  "reference_repeat_purity",
  "gencode_gene_region",
  "gencode_gene_name",
  "gencode_gene_id",
  "gencode_transcript_id",
  "refseq_gene_region",
  "mane_gene_region"
)

trexplorer_catalogue <- read_tsv(
  bed_file,
  col_types = cols(
    Chromosome = col_character(),
    start = col_integer(),
    end = col_integer(),
    motif = col_character(),
    motif_length = col_integer(),
    tract_length = col_integer(),
    copies = col_double(),
    locus_id = col_character(),
    variant_type = col_character(),
    source = col_character(),
    canonical_motif = col_character(),
    reference_repeat_purity = col_double(),
    gencode_gene_region = col_character(),
    gencode_gene_name = col_character(),
    gencode_gene_id = col_character(),
    gencode_transcript_id = col_character(),
    refseq_gene_region = col_character(),
    mane_gene_region = col_character()
  ),
  show_col_types = FALSE,
  progress = TRUE
)

trexplorer_filtered_1 <- trexplorer_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

trexplorer_removed_1 <- trexplorer_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(trexplorer_filtered_1, filtered_1_output_file)
write_tsv(trexplorer_removed_1, removed_1_output_file)

trexplorer_step_2 <- trexplorer_filtered_1 %>%
  mutate(
    minimum_tract_length = case_when(
      motif_length == 1L ~ 9L,
      motif_length == 2L ~ 10L,
      motif_length == 3L ~ 9L,
      motif_length == 4L ~ 12L,
      motif_length == 5L ~ 15L,
      motif_length == 6L ~ 18L,
      TRUE ~ NA_integer_
    )
  )

trexplorer_filtered_2 <- trexplorer_step_2 %>%
  filter(!is.na(tract_length), tract_length >= minimum_tract_length) %>%
  select(-minimum_tract_length)

trexplorer_removed_2 <- trexplorer_step_2 %>%
  filter(is.na(tract_length) | is.na(minimum_tract_length) | tract_length < minimum_tract_length) %>%
  select(-minimum_tract_length)

write_tsv(trexplorer_filtered_2, filtered_2_output_file)
write_tsv(trexplorer_removed_2, removed_2_output_file)