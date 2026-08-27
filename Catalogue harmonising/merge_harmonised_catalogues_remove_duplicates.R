library(tidyverse)

# duplicate = same chromosome, same motif length, and start coordinate within
#             +/- 1 bp of a previously encountered record.
# end coordinate is not used for duplicate status.
# only motif length is used.
# Records with  exact same start are same_start_duplicate_removed.
# Records with start coordinate +/- 1 bp are near_start_duplicate_removed.

base_dir <- Sys.getenv(
  "STR_CATALOGUE_BASE_DIR",
  unset = "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised"
)

input_files <- list.files(
  base_dir,
  pattern = "_harmonised\\.tsv$",
  full.names = TRUE
)

merged_output_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_with_duplicate_status_feature_preserving_v2.tsv"
)
deduplicated_output_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_deduplicated_feature_preserving_v2.tsv"
)
removed_duplicates_output_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_removed_duplicates_feature_preserving_v2.tsv"
)
duplicate_group_summary_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_duplicate_group_summary_feature_preserving_v2.tsv"
)
duplicate_counts_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_duplicate_counts_by_catalogue_feature_preserving_v2.tsv"
)
duplicate_counts_wide_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_duplicate_counts_by_catalogue_wide_feature_preserving_v2.tsv"
)
overall_summary_file <- file.path(
  base_dir,
  "merged_harmonised_catalogues_duplicate_summary_overall_feature_preserving_v2.tsv"
)

start_coordinate_tolerance_bp <- 1


motif %>%
    str_to_upper() %>%
    str_squish() # converts to upper case



read_harmonised_catalogue <- function(file_path) {
  read_tsv(
    file_path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>%
    mutate(source_file = basename(file_path), .before = 1L) # adds new column to say source file name
}


required_columns <- c(
  "catalogue", "chromosome", "start", "end", "motif",
  "motif_length", "repeat_length"
) # required columns 

harmonised_catalogues <- input_files %>%
  map_dfr(read_harmonised_catalogue) # reads all input catalogues and combines them into one table

missing_required_columns <- setdiff(required_columns, names(harmonised_catalogues))

if (length(missing_required_columns) > 0) {
  stop(
    "The harmonised files are missing required column(s): ",
    paste(missing_required_columns, collapse = ", ")
  ) # if missing required columns, code stop so I can check
}


merged_catalogue <- harmonised_catalogues %>%
  mutate(
    row_id = row_number(),
    chromosome = as.character(chromosome),
    start = parse_number(as.character(start)),
    end = parse_number(as.character(end)),
    motif = str_to_upper(as.character(motif)),
    motif_for_duplicate_check = map_chr(motif, clean_duplicate_motif),
    motif_length = coalesce(
      parse_number(as.character(motif_length)),
      if_else(
        !is.na(motif_for_duplicate_check),
        as.numeric(str_length(motif_for_duplicate_check)),
        NA_real_
      )
    ),
    repeat_length = parse_number(as.character(repeat_length)),
    locus_length = case_when(
      !is.na(repeat_length) & repeat_length > 0 ~ repeat_length,
      !is.na(start) & !is.na(end) ~ abs(end - start) + 1,
      TRUE ~ NA_real_
    )
  ) %>%
  relocate(
    row_id,
    source_file,
    catalogue,
    chromosome,
    start,
    end,
    motif,
    motif_for_duplicate_check,
    motif_length,
    repeat_length,
    locus_length
  ) #overall makes merged catalogue dataframe to clean and standardise columns. 
# ensures that some columns are numerical, capitals. then calculates motif length by counting characters in motif column
# then moves important columns to the beginning of the df

assign_duplicate_status <- function(group_data, group_key) {
  group_data <- group_data %>%
    arrange(start, catalogue, row_id)

  n_loci <- nrow(group_data)

  duplicate_status <- rep("unique", n_loci)
  removed_as_duplicate <- rep(FALSE, n_loci)
  local_duplicate_group <- rep(NA_integer_, n_loci)
  duplicate_of_row_id <- rep(NA_integer_, n_loci)
  duplicate_match_row_id <- rep(NA_integer_, n_loci)
  duplicate_overlap_bp <- rep(NA_real_, n_loci)
  duplicate_reciprocal_overlap <- rep(NA_real_, n_loci)
  duplicate_start_difference_bp <- rep(NA_real_, n_loci)
  duplicate_end_difference_bp <- rep(NA_real_, n_loci)
  possible_duplicate_of_row_id <- rep(NA_integer_, n_loci)
  possible_duplicate_reciprocal_overlap <- rep(NA_real_, n_loci)
  genomic_feature_original <- as.character(group_data$genomic_feature)
  genomic_feature_resolved <- genomic_feature_original
  genomic_feature_source_row_id <- if_else(
    !is_missing_genomic_feature(genomic_feature_original),
    group_data$row_id,
    NA_integer_
  )
  genomic_feature_inherited_from_duplicate <- rep(FALSE, n_loci)

  active_loci <- integer()
  group_representative <- integer()
  group_first_feature <- character()
  group_first_feature_source_index <- integer()
  next_group_number <- 0L

  for (i in seq_len(n_loci)) {
    if (length(active_loci) > 0) {
      active_loci <- active_loci[
        group_data$start[active_loci] >=
          group_data$start[[i]] - start_coordinate_tolerance_bp
      ]
    }

    best_duplicate_index <- NA_integer_
    best_start_difference <- NA_real_

    if (length(active_loci) > 0) {
      start_difference <- abs(group_data$start[active_loci] - group_data$start[[i]])
      duplicate_hits <- which(start_difference <= start_coordinate_tolerance_bp)

      if (length(duplicate_hits) > 0) {
        best_duplicate_hit <- duplicate_hits[[which.min(start_difference[duplicate_hits])]]
        best_duplicate_index <- active_loci[[best_duplicate_hit]]
        best_start_difference <- start_difference[[best_duplicate_hit]]
      }
    }

    if (!is.na(best_duplicate_index)) {
      previous_group <- local_duplicate_group[[best_duplicate_index]]

      local_duplicate_group[[i]] <- previous_group
      representative_index <- group_representative[[previous_group]]

      duplicate_status[[i]] <- if_else(
        best_start_difference == 0,
        "same_start_duplicate_removed",
        "near_start_duplicate_removed"
      )
      removed_as_duplicate[[i]] <- TRUE
      duplicate_of_row_id[[i]] <- group_data$row_id[[representative_index]]
      duplicate_match_row_id[[i]] <- group_data$row_id[[best_duplicate_index]]
      duplicate_start_difference_bp[[i]] <- best_start_difference
    } else {
      next_group_number <- next_group_number + 1L
      local_duplicate_group[[i]] <- next_group_number
      group_representative[[next_group_number]] <- i
    }

    current_group <- local_duplicate_group[[i]]

    if (
      !is_missing_genomic_feature(genomic_feature_original[[i]]) &&
      (
        length(group_first_feature) < current_group ||
          is.na(group_first_feature[[current_group]])
      )
    ) {
      group_first_feature[[current_group]] <- genomic_feature_original[[i]]
      group_first_feature_source_index[[current_group]] <- i
    }

    active_loci <- c(active_loci, i)
  }

  for (group_number in seq_len(next_group_number)) {
    representative_index <- group_representative[[group_number]]

    if (
      is_missing_genomic_feature(
        genomic_feature_original[[representative_index]]
      ) &&
      length(group_first_feature) >= group_number &&
      !is.na(group_first_feature[[group_number]])
    ) {
      donor_index <- group_first_feature_source_index[[group_number]]
      genomic_feature_resolved[[representative_index]] <-
        genomic_feature_original[[donor_index]]
      genomic_feature_source_row_id[[representative_index]] <-
        group_data$row_id[[donor_index]]
      genomic_feature_inherited_from_duplicate[[representative_index]] <-
        donor_index != representative_index
    }
  }

  group_data %>%
    mutate(
      duplicate_group_id = str_c(
        group_key$chromosome[[1]],
        group_key$motif_length[[1]],
        local_duplicate_group,
        sep = "|"
      ),
      duplicate_status = duplicate_status,
      removed_as_duplicate = removed_as_duplicate,
      duplicate_of_row_id = duplicate_of_row_id,
      duplicate_match_row_id = duplicate_match_row_id,
      duplicate_overlap_bp = duplicate_overlap_bp,
      duplicate_reciprocal_overlap = duplicate_reciprocal_overlap,
      duplicate_start_difference_bp = duplicate_start_difference_bp,
      duplicate_end_difference_bp = duplicate_end_difference_bp,
      possible_duplicate_of_row_id = possible_duplicate_of_row_id,
      possible_duplicate_reciprocal_overlap = possible_duplicate_reciprocal_overlap,
      genomic_feature_original = genomic_feature_original,
      genomic_feature_resolved = genomic_feature_resolved,
      genomic_feature_source_row_id = genomic_feature_source_row_id,
      genomic_feature_inherited_from_duplicate =
        genomic_feature_inherited_from_duplicate
    )
} # v large chunk to compare repeat w/ chromosome + motif length groups once sorting by start coord
#transitive grouping


catalogue_ready_for_duplicate_check <- merged_catalogue %>%
  filter(
    !is.na(chromosome),
    !is.na(start),
    !is.na(motif_length)
  )

catalogue_with_duplicate_status <- catalogue_ready_for_duplicate_check %>%
  group_by(chromosome, motif_length) %>%
  group_modify(assign_duplicate_status) %>%
  ungroup() %>%
  bind_rows(catalogue_not_checked) %>%
  arrange(row_id) #divides by chromosome + motif length

deduplicated_catalogue <- catalogue_with_duplicate_status %>%
  filter(!removed_as_duplicate) %>%
  mutate(genomic_feature = genomic_feature_resolved) %>%
  select(-genomic_feature_resolved) # remove duplicate repeats, jeep unique repeats

removed_duplicates <- catalogue_with_duplicate_status %>%
  filter(removed_as_duplicate) # in the name = stores the duplicates to a catalogue



duplicate_group_summary <- catalogue_with_duplicate_status %>%
  filter(!is.na(duplicate_group_id)) %>%
  group_by(duplicate_group_id, chromosome, motif_length) %>%
  summarise(
    n_records = n(),
    n_removed_duplicates = sum(removed_as_duplicate),
    n_retained_records = sum(!removed_as_duplicate),
    n_catalogues = n_distinct(catalogue),
    catalogues_present = str_c(sort(unique(catalogue)), collapse = "; "),
    motifs_present = str_c(sort(unique(na.omit(motif))), collapse = "; "),
    min_start = min(start, na.rm = TRUE),
    max_start = max(start, na.rm = TRUE),
    start_span_bp = max(start, na.rm = TRUE) - min(start, na.rm = TRUE),
    max_end = max(end, na.rm = TRUE),
    representative_row_id = row_id[which(!removed_as_duplicate)[[1]]],
    duplicate_statuses = str_c(sort(unique(duplicate_status)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_records > 1 | n_removed_duplicates > 0)
#table of how many and which repeats wer in a duplicate group


duplicate_counts_by_catalogue <- catalogue_with_duplicate_status %>%
  count(catalogue, duplicate_status, removed_as_duplicate, name = "n_loci") %>%
  arrange(catalogue, duplicate_status)
# how many repeats were removed per catlogue


duplicate_counts_wide <- duplicate_counts_by_catalogue %>%
  select(catalogue, duplicate_status, n_loci) %>%
  pivot_wider(
    names_from = duplicate_status,
    values_from = n_loci,
    values_fill = 0
  )
#makes table to clean up previous part


overall_summary <- tibble(
  metric = c(
    "input_harmonised_files",
    "merged_rows",
    "deduplicated_rows",
    "removed_duplicate_rows",
    "same_start_duplicate_removed_rows",
    "near_start_duplicate_removed_rows",
    "deduplicated_rows_inheriting_genomic_feature",
    "not_checked_missing_required_value_rows",
    "start_coordinate_tolerance_bp",
    "deduplication_chromosome_rule",
    "deduplication_motif_rule",
    "deduplication_end_coordinate_rule"
  ),
  value = as.character(c(
    length(input_files),
    nrow(catalogue_with_duplicate_status),
    nrow(deduplicated_catalogue),
    nrow(removed_duplicates),
    sum(catalogue_with_duplicate_status$duplicate_status == "same_start_duplicate_removed"),
    sum(catalogue_with_duplicate_status$duplicate_status == "near_start_duplicate_removed"),
    sum(
      deduplicated_catalogue$genomic_feature_inherited_from_duplicate,
      na.rm = TRUE
      ),
    )
  )
)
#summary of input and output counts. categories of duplicate





# Write to new v2 filenames. Existing outputs are left untouched.
write_tsv(deduplicated_catalogue, deduplicated_output_file)
verify_written_rows(deduplicated_output_file, nrow(deduplicated_catalogue))

write_tsv(catalogue_with_duplicate_status, merged_output_file)
verify_written_rows(
  merged_output_file,
  nrow(catalogue_with_duplicate_status)
)

write_tsv(removed_duplicates, removed_duplicates_output_file)
verify_written_rows(
  removed_duplicates_output_file,
  nrow(removed_duplicates)
)

write_tsv(duplicate_group_summary, duplicate_group_summary_file)
write_tsv(duplicate_counts_by_catalogue, duplicate_counts_file)
write_tsv(duplicate_counts_wide, duplicate_counts_wide_file)
write_tsv(overall_summary, overall_summary_file)
