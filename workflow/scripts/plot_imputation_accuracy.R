suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))


custom_colors <- c("gray", "black", "#1F77B4", "#9467BD", "#2CA02C", "#D62728", "#FF7F0E")
names(custom_colors) <- c("LCS", "TRUTH", "BEAGLE", "GENEIMP", "GLIMPSE", "QUILT", "STITCH")

# KGD
# kgd_lcs_file <- snakemake@input[["KGD_LCS"]]
# kgd_imputed_files <- snakemake@input[["KGD_IMPUTED"]]
# kgd_truth_file <- snakemake@input[["KGD_TRUTH"]]

# kgd_lcs_file <- "/home/tvi069/imputation/test/OUTPUT/KGD/LCS/p4_G5.txt"
# kgd_truth_files<- "/home/tvi069/imputation/test/OUTPUT/KGD/TRUTH/p4_G5.txt"
# kgd_imputed_file <- list.files("/home/tvi069/imputation/test/OUTPUT/KGD/IMPUTED", 
#                                "G5.txt", recursive = T, full.names = T)

# kgd_lcs <- data.frame(G5 = as.numeric(readLines(kgd_lcs_file)),
#                       data = "LCS")
# kgd_truth <- data.frame(G5 = as.numeric(readLines(kgd_truth_file)),
#                       data = "TRUTH")
# kgd_imputed <- do.call(rbind, lapply(kgd_imputed_files, function(f) {
#   tool <- basename(dirname(f))
#   data.frame(G5 = as.numeric(readLines(f)),
#                     data = tool)
# }))

# kgd <- do.call(rbind, list(kgd_lcs, kgd_truth, kgd_imputed))
# kgd$data <- factor(kgd$data, c("LCS", "TRUTH", sort(unique(kgd_imputed$data))))

# png(snakemake@output[["KGD"]], width = 1500, height = 500)
# ggplot(kgd) +
#   facet_grid(cols = vars(data)) +
#   geom_histogram(aes(x = G5, fill = data), bins = 20, show.legend = F) + 
#   theme_minimal() + 
#   theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         strip.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12)) +
#   scale_color_manual(values = custom_colors) +
#   scale_fill_manual(values = custom_colors)
# dev.off()


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
nr$data <- factor(nr$data, c("LCS", "TRUTH", sort(unique(nr_imputed$data))))

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
  xlab("Pairwise Relatedness (rab)") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()

png(snakemake@output[["NGSRELATE_INBREEDING"]], width = 1500, height = 500)
nr %>% 
  dplyr::group_by(a, data) %>% 
  dplyr::summarise(mean_Fa = mean(Fa)) %>% 
  ggplot() +
  facet_grid(cols =vars(data)) +
  geom_histogram(aes(x = mean_Fa, fill = data), bins = 20, show.legend = F) + 
  theme_minimal() + 
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  xlab("Individual Inbreeding Coefficient (Fa)") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()


# # IQS
# iqs_var_files <- snakemake@input[["IQS_VAR"]]
# # iqs_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "_per_variant_results.txt",
# #                         recursive = T, full.names = T)

# iqs_var <- do.call(rbind, lapply(iqs_var_files, function(f) {
#   tool <- basename(dirname(f))
#   i <- read.delim(f)
#   i$data <- tool
#   i
# }))

# iqs_var <- iqs_var %>% 
#   pivot_longer(cols = c("IQS", "precision", "recall", "F.score"),
#                names_to = "Metrics", values_to = "Value") 
# iqs_var$MAF_bin <- cut(iqs_var$WGS_MAF, breaks = c(0, 0.05, 0.1, 0.25, 0.5), include.lowest = T)
# iqs_var$Metrics <- factor(iqs_var$Metrics, levels = c("recall", "precision", "F.score", "IQS"))

# png(snakemake@output[["IQS_VAR"]], width = 1500, height = 800)
# ggplot(iqs_var) +
#   facet_grid(cols = vars(data), rows = vars(Metrics)) +
#   geom_boxplot(aes(y = Value, x = MAF_bin, color = data), outliers = F, show.legend = F) +
#   theme_minimal() +
#   theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         strip.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12)) + 
#   ylab("Value per variant") + 
#   xlab("Truth MAF") +
#   scale_color_manual(values = custom_colors) +
#   scale_fill_manual(values = custom_colors)
# dev.off()


# iqs_sam_files <- snakemake@input[["IQS_SAM"]]
# # iqs_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "_per_variant_results.txt",
# #                         recursive = T, full.names = T)

# iqs_sam <- do.call(rbind, lapply(iqs_sam_files, function(f) {
#   tool <- basename(dirname(f))
#   i <- read.delim(f)
#   i$data <- tool
#   i
# }))

# iqs_sam <- iqs_sam %>% 
#   pivot_longer(cols = c("precision", "recall", "F.score"),
#                names_to = "Metrics", values_to = "Value") 
# iqs_sam$Metrics <- factor(iqs_sam$Metrics, levels = c("recall", "precision", "F.score"))

# png(snakemake@output[["IQS_SAM"]], width = 800, height = 800)
# ggplot(iqs_sam) +
#   facet_grid(rows = vars(Metrics)) +
#   geom_boxplot(aes(y = Value, x = data, color = data), outliers = F, show.legend = F) +
#   theme_minimal() +
#   theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
#         strip.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12)) + 
#   ylab("Value per sample") + 
#   xlab("Imputation") +
#   scale_color_manual(values = custom_colors) +
#   scale_fill_manual(values = custom_colors)
# dev.off()


# Hap.py
happy_files <- snakemake@input[["HAPPY"]]
# happy_files <- list.files("/nesi/nobackup/uoa04053/tram/hihi/imputation_2autosomes/OUTPUT/ACCURACY",
#                           ".+_[A-Z]+_accuracy.summary.csv", recursive = T, full.names = T)


happy <- do.call(rbind, lapply(happy_files, function(f) {
  # print(f)
  hp <- read.csv(f)[1,]
  tool <- basename(dirname(f))
  hp$data <- tool
  hp
}))

happy <- happy %>% 
  pivot_longer(cols = c("METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"),
               names_to = "Metrics", values_to = "Value") 
happy$Metrics <- factor(happy$Metrics, levels = c("METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"))

png(snakemake@output[["HAPPY"]], width = 800, height = 800)
ggplot(happy) +
  facet_grid(rows = vars(Metrics), scales = "free_x") +
  geom_col(aes(y = Value, x = data, fill = data), show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylab("Value") + 
  xlab("Imputation") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()

# Customized accuracy metrics
var_files <- snakemake@input[["METRICS_VAR"]]
# var_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "accuracy_metrics_per_variant.tsv",
#                         recursive = T, full.names = T)

var <- do.call(rbind, lapply(var_files, function(f) {
  tool <- basename(dirname(f))
  v <- read.delim(f)
  v$data <- tool
  v
}))


var <- var %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQS"),
               names_to = "Metrics", values_to = "Value") 
var$Truth_MAF <- sapply(var$Truth_ALT_Freq, function(x) min(c(x, 1-x)))
var$MAF_bin <- cut(var$Truth_MAF, breaks = c(0, 0.05, 0.1, 0.25, 0.5), include.lowest = T)
var$Metrics <- factor(var$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Recall", "Precision", "F1", "IQS"))

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
  xlab("Truth MAF") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()


sam_files <- snakemake@input[["METRICS_SAM"]]
# sam_files <- list.files("/home/tvi069/imputation/test/OUTPUT/ACCURACY", "accuracy_metrics_per_sample.tsv",
#                         recursive = T, full.names = T)

sam <- do.call(rbind, lapply(sam_files, function(f) {
  tool <- basename(dirname(f))
  s <- read.delim(f)
  s$data <- tool
  s
}))

sam <- sam %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQS"),
               names_to = "Metrics", values_to = "Value") 
sam$Metrics <- factor(sam$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQ"))

png(snakemake@output[["SAMPLE"]], width = 1000, height = 800)
ggplot(sam) +
  facet_grid(rows = vars(Metrics), scales = "free_x") +
  geom_boxplot(aes(y = Value, x = data, color = data), outliers = F, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  ylab("Value per sample") + 
  xlab("Imputation") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()

# MAF
var$Imputed_MAF <- sapply(var$Imputed_ALT_Freq, function(x) min(c(x, 1-x)))

png(snakemake@output[["MAF"]], width = 1500, height = 300)
ggplot(var, aes(y = Imputed_MAF, x = Truth_MAF, color = data)) +
  facet_grid(cols = vars(data)) +
  geom_point(alpha = 0.7, show.legend = F) +
  geom_smooth(method = lm, se = T, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) + 
  xlab("Groundtruth MAF") + 
  ylab("Imputed MAF") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()


# Benchmark

benchmark_files <- list.files(snakemake@params[["BENCHMARK"]], full.names = T)
# benchmark_files <- list.files("/home/tvi069/imputation/test/OUTPUT/BENCHMARK",
#                               full.names = T)

benchmark <- do.call(rbind, lapply(benchmark_files, function(f) {
  b <- data.frame(Imputation = gsub(".txt", "", basename(f)),
                  read.delim(f))
  # b$Imputation <- gsub(".txt", "", basename(f))
  b
}))

png(snakemake@output[["BENCHMARK"]], width=1000, height=300)
p<-tableGrob(benchmark)
grid.arrange(p)
dev.off()
