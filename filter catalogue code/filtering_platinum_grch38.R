library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/PlatinumTRs"
bed_file <- file.path(base_dir, "human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed")

filtered_1_output_file <- file.path(base_dir, "platinum_grch38_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "platinum_grch38_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "platinum_grch38_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "platinum_grch38_removed_2.0.tsv")

catalogue_name <- "platinum_palladium_GRCh38"

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)

platinum_catalogue <- read_tsv(
  bed_file,
  quote = "",
  col_types = cols(
    `#chrom` = col_character(),
    start = col_integer(),
    end = col_integer(),
    TRid = col_character(),
    longest_homopolymer = col_integer(),
    gc_content = col_double(),
    n_motifs = col_integer(),
    min_motiflen = col_integer(),
    max_motiflen = col_integer(),
    motifs = col_character(),
    struc = col_character(),
    SegDup = col_integer(),
    Telomere = col_integer(),
    Centromere = col_integer(),
    CDS = col_integer(),
    .default = col_skip()
  )
) %>%
  transmute(
    catalogue = catalogue_name,
    tr_id = TRid,
    longest_homopolymer = longest_homopolymer,
    gc_content = gc_content,
    n_motifs = n_motifs,
    min_motif_length = min_motiflen,
    max_motif_length = max_motiflen,
    region_chrom = `#chrom`,
    region_start = start,
    region_end = end,
    motif = str_split(motifs, ","),
    locus_structure = struc,
    segmental_duplication_overlap = SegDup,
    telomere_overlap = Telomere,
    centromere_overlap = Centromere,
    cds_overlap = CDS
  ) %>%
  unnest(motif) %>%
  mutate(
    motif_length = nchar(motif),
    repeat_length = region_end - region_start,
    motif_copies = repeat_length / motif_length
  ) %>%
  select(
    catalogue,
    tr_id,
    longest_homopolymer,
    gc_content,
    n_motifs,
    min_motif_length,
    max_motif_length,
    region_chrom,
    region_start,
    region_end,
    motif,
    motif_length,
    motif_copies,
    repeat_length,
    locus_structure,
    segmental_duplication_overlap,
    telomere_overlap,
    centromere_overlap,
    cds_overlap
)

platinum_filtered_1 <- platinum_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

platinum_removed_1 <- platinum_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(platinum_filtered_1, filtered_1_output_file)
write_tsv(platinum_removed_1, removed_1_output_file)


platinum_step_2 <- platinum_filtered_1 %>%
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

platinum_filtered_2 <- platinum_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

platinum_removed_2 <- platinum_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(platinum_filtered_2, filtered_2_output_file)
write_tsv(platinum_removed_2, removed_2_output_file)

