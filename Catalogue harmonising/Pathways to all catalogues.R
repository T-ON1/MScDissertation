library(tidyverse)

illumina <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/illumina 174K/illumina_filtered_2.0.tsv")
UCSC <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/UCSC/ucsc_filtered_2.0.tsv")
VAMOS <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/VAMOS original VCF/vamos_filtered_2.0.tsv")
gangSTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/GangSTR-master/gangstr_filtered_2.0.tsv")
hipSTR <-  read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/HipSTR/hipstr_filtered_2.0.tsv")
adotto <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/adottoTR_filtering/adotto_TRregions_v1.2_filtered_2.0.tsv")
popSTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/popSTR 5,401,401/popstr_filtered_2.0.tsv")
platinumTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/PlatinumTRs/platinum_grch38_filtered_2.0.tsv")
chiu <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Chiu catalog/chiu_hg38.v1_filtered_2.0.tsv")
TRExplorer <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/TRexplorer/trexplorer_filtered_2.0.tsv")


#harmonising names of columns
# reordering columns to the same format
#catalogue, chromosome, start, end, motif, motif_length, motif_copies, repeat_length

#removed removed nothing, renamed headers
harmonised_illumina <- illumina %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         calculated_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length,
         reference_region,
         locus_structure,
         locus_id,
         variant_id,
         variant_type
)

#removed ucsc_name and one of the motif copies, renamed headers
harmonised_UCSC <- UCSC %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         source_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         source_motif_copies,
         calculated_motif_copies,
         repeat_length
)

#removed vcf_id, vcf_filter, svtype, renamed headers
harmonised_VAMOS <- VAMOS %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         calculated_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length
)


#removed repeat_sequence, renamed headers
harmonised_GangSTR <- gangSTR %>% 
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         calculated_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length
)

#removed locus_id it was just "Human_STR_Numbers", it wasnt a gene. removed one of the motif lengths, renamed headers
harmonised_HipSTR <- hipSTR %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         source_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         source_motif_copies,
         calculated_motif_copies,
         repeat_length
)


#renamed columns and rearranged them 
harmonised_adotto <- adotto %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         source_motif_copies = motif_copies,
         variant_type = repeat_class)%>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         source_motif_copies,
         calculated_motif_copies,
         repeat_length,
         genomic_feature,
         variant_type
)

  
#removed irrelevant columns e.g. slippage, stutter, fraction of ATCG. renamed headers
harmonised_PopSTR <- popSTR %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         source_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         source_motif_copies,
         calculated_motif_copies,
         repeat_length
)


#removed irrelevant columns e.g. gc content + overlap, renamed + rearranged headers  
harmonised_PlatinumTR <- platinumTR %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         calculated_motif_copies = motif_copies) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length,
         locus_structure
)

  
#rearranged + renamed headers
harmonised_chiu <- chiu %>%
  rename(chromosome = region_chrom,
         start = region_start,
         end = region_end,
         calculated_motif_copies = motif_copies,
         genomic_feature = feature) %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length,
         genomic_feature
)

  
#added catalogue column, renamed columns, + rearranged headers 
harmonised_TRExplorer <- TRExplorer %>%
  rename(chromosome = Chromosome,
         repeat_length = tract_length,
         calculated_motif_copies = copies) %>%
  mutate(catalogue = "TRExplorer") %>%
  select(catalogue,
         chromosome,
         start,
         end,
         motif,
         motif_length,
         calculated_motif_copies,
         repeat_length
)


#writing all of the harmonised dataframes to the harmonised folder

write_tsv(harmonised_illumina, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/illumina_harmonised.tsv")
write_tsv(harmonised_UCSC, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/UCSC_harmonised.tsv")
write_tsv(harmonised_VAMOS, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/VAMOS_harmonised.tsv")
write_tsv(harmonised_GangSTR, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/GangSTR_harmonised.tsv")
write_tsv(harmonised_HipSTR, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/HipSTR_harmonised.tsv")
write_tsv(harmonised_adotto, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/adotto_harmonised.tsv")
write_tsv(harmonised_PopSTR, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/PopSTR_harmonised.tsv")
write_tsv(harmonised_PlatinumTR, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/PlatinumTR_harmonised.tsv")
write_tsv(harmonised_chiu, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/chiu_harmonised.tsv")
write_tsv(harmonised_TRExplorer, "C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/TRExplorer_harmonised.tsv")


illumina <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/illumina_harmonised.tsv")
UCSC <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/UCSC_harmonised.tsv")
VAMOS <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/VAMOS_harmonised.tsv")
gangSTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/GangSTR_harmonised.tsv")
hipSTR <-  read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/HipSTR_harmonised.tsv")
adotto <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/adotto_harmonised.tsv")
popSTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/PopSTR_harmonised.tsv")
platinumTR <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/PlatinumTR_harmonised.tsv")
chiu <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/chiu_harmonised.tsv")
TRExplorer <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/Catalogues_harmonised/TRExplorer_harmonised.tsv")













plat <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/PlatinumTRs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed")

plat_raw <- read_tsv("C:/Users/tiern/Desktop/Dissertation catalogues/PlatinumTRs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations_filtered.tsv")


plat_composite <- plat_raw %>%
  filter(grepl(")n(", a, fixed = TRUE))









