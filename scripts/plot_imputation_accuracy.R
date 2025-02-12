suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# KGD
kgd_lcs_file <- snakemake@input[["KGD_LCS"]]
kgd_imputed_files <- snakemake@input[["KGD_IMPUTED"]]
kgd_truth_file <- snakemake@input[["KGD_TRUTH"]]

# kgd_lcs_file <- "/home/tvi069/imputation/test/OUTPUT/KGD/LCS/p4_G5.txt"
# kgd_truth_files<- "/home/tvi069/imputation/test/OUTPUT/KGD/TRUTH/p4_G5.txt"
# kgd_imputed_file <- list.files("/home/tvi069/imputation/test/OUTPUT/KGD/IMPUTED", 
#                                "G5.txt", recursive = T, full.names = T)

kgd_lcs <- data.frame(G5 = as.numeric(readLines(kgd_lcs_file)),
                      data = "LCS")
kgd_truth <- data.frame(G5 = as.numeric(readLines(kgd_truth_file)),
                      data = "TRUTH")
kgd_imputed <- do.call(rbind, lapply(kgd_imputed_files, function(f) {
  tool <- basename(dirname(f))
  data.frame(G5 = as.numeric(readLines(f)),
                    data = tool)
}))

kgd <- do.call(rbind, list(kgd_lcs, kgd_truth, kgd_imputed))
kgd$data <- factor(kgd$data, c("LCS", "TRUTH", unique(kgd_imputed$data)))

png(snakemake@output[["KGD"]], width = 1500, height = 500)
ggplot(kgd) +
  facet_grid(cols = vars(data)) +
  geom_histogram(aes(x = G5, fill = data), bins = 20, show.legend = F) + 
  theme_minimal() + 
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
dev.off()


# NGSRELATE
nr_lcs_file <- snakemake@input[["NGSRELATE_BAM"]]
nr_imputed_files <- snakemake@input[["NGSRELATE_IMPUTED"]]
nr_truth_file <- snakemake@input[["NGSRELATE_TRUTH"]]

# nr_lcs_file <- "/home/tvi069/imputation/test/OUTPUT/NGSRELATE/BAM/p4.ngsrelate"
# nr_imputed_files <- list.files("/home/tvi069/imputation/test/OUTPUT/NGSRELATE/IMPUTED", 
#                                ".ngsrelate", recursive = T, full.names = T)
# nr_truth_file <- "/home/tvi069/imputation/test/OUTPUT/NGSRELATE/TRUTH/p4.ngsrelate"

nr_lcs <- read.delim(nr_lcs_file)
nr_lcs$data <- "LCS"
nr_truth <- read.delim(nr_truth_file)
nr_truth$data <- "TRUTH"
nr_imputed <- do.call(rbind, lapply(nr_imputed_files, function(f) {
  tool <- gsub(".+_|\\.ngsrelate", "", basename(f))
  nr <- read.delim(f)
  nr$data <- tool
  nr
}))

nr <- do.call(rbind, list(nr_lcs, nr_truth, nr_imputed))
nr$data <- factor(nr$data, c("LCS", "TRUTH", unique(nr_imputed$data)))

png(snakemake@output[["NGSRELATE_REL"]], width = 1500, height = 500)
ggplot(nr) +
  facet_grid(cols = vars(data)) +
  geom_histogram(aes(x = rab, fill = data), bins = 20, show.legend = F) + 
  theme_minimal() + 
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  xlab("Pairwise Relatedness (rab)")
dev.off()

png(snakemake@output[["NGSRELATE_INBREEDING"]], width = 1500, height = 500)
nr %>% 
  dplyr::group_by(a, data) %>% 
  dplyr::summarise(mean_Fa = mean(Fa)) %>% 
  ggplot() +
  facet_grid(cols = vars(data)) +
  geom_histogram(aes(x = mean_Fa, fill = data), bins = 20, show.legend = F) + 
  theme_minimal() + 
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  xlab("Individual Inbreeding Coefficient (Fa)")
dev.off()


# IQS
iqs_files <- snakemake@input[["IQS_VAR"]]
# iqs_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "_per_variant_results.txt",
#                         recursive = T, full.names = T)

iqs <- do.call(rbind, lapply(iqs_files, function(f) {
  tool <- basename(dirname(f))
  i <- read.delim(f)
  i$data <- tool
  i
}))

iqs <- iqs %>% 
  pivot_longer(cols = c("IQS", "precision", "recall", "F.score"),
               names_to = "Metrics", values_to = "Value") 
iqs$MAF_bin <- cut(iqs$WGS_MAF, breaks = c(0, 0.05, 0.1, 0.25, 0.5), include.lowest = T)
iqs$Metrics <- factor(iqs$Metrics, levels = c("recall", "precision", "F.score", "IQS"))

png(snakemake@output[["IQS"]], width = 1500, height = 800)
ggplot(iqs) +
  facet_grid(cols = vars(data), rows = vars(Metrics)) +
  geom_boxplot(aes(y = Value, x = MAF_bin, color = data), outliers = F, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylab("Value per variant") + 
  xlab("Truth MAF")
dev.off()

# Hap.py


# Customized accuracy metrics
var_files <- snakemake@input[["METRICS_VAR"]]
var_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "accuracy_metrics_per_variant.tsv",
                        recursive = T, full.names = T)

var <- do.call(rbind, lapply(var_files, function(f) {
  tool <- basename(dirname(f))
  v <- read.delim(f)
  v$data <- tool
  v
}))


var <- var %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1"),
               names_to = "Metrics", values_to = "Value") 
var$Truth_MAF <- sapply(var$Truth_ALT_Freq, function(x) min(c(x, 1-x)))
var$MAF_bin <- cut(var$Truth_MAF, breaks = c(0, 0.05, 0.1, 0.25, 0.5), include.lowest = T)
var$Metrics <- factor(var$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1"))

png(snakemake@output[["VARIANT"]], width = 1600, height = 900)
ggplot(var) +
  facet_grid(cols = vars(data), rows = vars(Metrics)) +
  geom_boxplot(aes(y = Value, x = MAF_bin, color = data), outliers = F, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylab("Value per variant") + 
  xlab("Truth MAF")
dev.off()


sam_files <- snakemake@input[["METRICS_SAM"]]
sam_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "accuracy_metrics_per_sample.tsv",
                        recursive = T, full.names = T)

sam <- do.call(rbind, lapply(sam_files, function(f) {
  tool <- basename(dirname(f))
  s <- read.delim(f)
  s$data <- tool
  s
}))

sam <- sam %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1"),
               names_to = "Metrics", values_to = "Value") 
sam$Metrics <- factor(sam$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1"))

png(snakemake@output[["SAMPLE"]], width = 1000, height = 800)
ggplot(sam) +
  facet_grid(cols = vars(data), rows = vars(Metrics), scales = "free_x") +
  geom_boxplot(aes(y = Value, x = data, color = data), outliers = F, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylab("Value per sample") + 
  xlab("Imputation")
dev.off()
