library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

accuracy_file <- args[1]
# accuracy_file <- "/nesi/nobackup/uoa04053/tram/sim_random/imputation/p4_sub20_rate080/50gen/target1235/ACCURACY/GLIMPSE/p4_sub20_rate080_target1235_GlimpseImputed_per_variant_results_test.txt"
wgs_acc <- fread(accuracy_file)
# wgs_acc <- wgs_acc[, -21]
wgs_acc <- wgs_acc %>% 
  mutate(CHROM = gsub(":\\d+", "", position),
         POS = gsub(".+:", "", position)) %>% 
  select(-position) %>% 
  relocate(CHROM, POS)

freq_file <- args[2]
# freq_file <- "/nesi/nobackup/uoa04053/tram/sim_random/imputation/p4_sub20_rate080/50gen/target1235/ACCURACY/BEAGLE/truth.frq"
freq <- fread(freq_file, skip = 1, header = F)
colnames(freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ")
freq <- freq %>% 
  mutate(CHROM = as.character(CHROM),
         POS = as.character(POS))

wgs_acc <- full_join(wgs_acc, freq, by = c("CHROM", "POS"))
# wgs_acc %>% filter(is.na(SNP))

fwrite(wgs_acc, accuracy_file, sep = "\t")
