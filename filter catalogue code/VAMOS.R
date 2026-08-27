library(tidyverse)

# Filter v1.0 keeps STR motifs with motif length 1-6 bp.
# Filter v2.0 starts from the v1.0 filtered output and applies motif-specific
# minimum repeat lengths:
# 1 bp motifs need >= 9 bp, 2 bp motifs need >= 10 bp, 3 bp motifs need >= 9 bp,
# 4 bp motifs need >= 12 bp, 5 bp motifs need >= 15 bp, and 6 bp motifs need >= 18 bp.

base_dir <- "C:/Users/tiern/Desktop/Dissertation catalogues/VAMOS original VCF"
vcf_file <- file.path(base_dir, "samples_rmDup_GRCh38_ori.vcf")

filtered_1_output_file <- file.path(base_dir, "vamos_filtered_1.0.tsv")
removed_1_output_file <- file.path(base_dir, "vamos_removed_1.0.tsv")
filtered_2_output_file <- file.path(base_dir, "vamos_filtered_2.0.tsv")
removed_2_output_file <- file.path(base_dir, "vamos_removed_2.0.tsv")

catalogue_name <- "vamos_GRCh38_ori"

step_1_output_files <- c(filtered_1_output_file, removed_1_output_file)
step_2_output_files <- c(filtered_2_output_file, removed_2_output_file)


vcf_header <- read_lines(vcf_file, n_max = 50)
skip_lines <- which(str_detect(vcf_header, "^#CHROM")) - 1

vamos_raw <- read_tsv(
  vcf_file,
  skip = skip_lines,
  col_types = cols_only(
    `#CHROM` = col_character(),
    POS = col_integer(),
    ID = col_character(),
    FILTER = col_character(),
    INFO = col_character()
  ),
  show_col_types = FALSE,
  progress = TRUE
)

vamos_catalogue <- vamos_raw %>%
  transmute(
    catalogue = catalogue_name,
    vcf_id = ID,
    vcf_filter = FILTER,
    region_chrom = `#CHROM`,
    region_start = POS,
    region_end = str_extract(INFO, "END=[^;]+"),
    motif = str_extract(INFO, "RU=[^;]+"),
    svtype = str_extract(INFO, "SVTYPE=[^;]+")
  ) %>%
  mutate(
    region_end = as.integer(str_remove(region_end, "END=")),
    motif = str_remove(motif, "RU="),
    motif = str_split(motif, ","),
    svtype = str_remove(svtype, "SVTYPE=")
  ) %>%
  unnest(motif) %>%
  mutate(
    motif_length = nchar(motif),
    repeat_length = region_end - region_start + 1,
    motif_copies = repeat_length / motif_length
  ) %>%
  select(
    catalogue,
    vcf_id,
    vcf_filter,
    region_chrom,
    region_start,
    region_end,
    motif,
    motif_length,
    motif_copies,
    repeat_length,
    svtype
  )

vamos_filtered_1 <- vamos_catalogue %>%
  filter(!is.na(motif_length), between(motif_length, 1L, 6L))

vamos_removed_1 <- vamos_catalogue %>%
  filter(is.na(motif_length) | !between(motif_length, 1L, 6L))

write_tsv(vamos_filtered_1, filtered_1_output_file)
write_tsv(vamos_removed_1, removed_1_output_file)


vamos_step_2 <- vamos_filtered_1 %>%
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

vamos_filtered_2 <- vamos_step_2 %>%
  filter(!is.na(repeat_length), repeat_length >= minimum_repeat_length) %>%
  select(-minimum_repeat_length)

vamos_removed_2 <- vamos_step_2 %>%
  filter(is.na(repeat_length) | is.na(minimum_repeat_length) | repeat_length < minimum_repeat_length) %>%
  select(-minimum_repeat_length)

write_tsv(vamos_filtered_2, filtered_2_output_file)
write_tsv(vamos_removed_2, removed_2_output_file)