suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))


custom_colors <- c("gray", "black", "#1F77B4", "#9467BD", "#2CA02C", "#D62728", "#FF7F0E")
names(custom_colors) <- c("LCS", "TRUTH", "BEAGLE", "GENEIMP", "GLIMPSE", "QUILT", "STITCH")


# NGSRELATE
nr_lcs_file <- snakemake@input[["NGSRELATE_BAM"]]
nr_imputed_file <- snakemake@input[["NGSRELATE_IMPUTED"]]
nr_truth_file <- snakemake@input[["NGSRELATE_TRUTH"]]

nr_lcs <- read.delim(nr_lcs_file)
nr_lcs$data <- "LCS"
nr_truth <- read.delim(nr_truth_file)
nr_truth$data <- "TRUTH"

nr_imputed <- read.delim(nr_imputed_file)
nr <- do.call(rbind, list(nr_lcs, nr_truth, nr_imputed))
nr$data <- factor(nr$data, c("LCS", "TRUTH", sort(unique(nr_imputed$data))))

png(snakemake@output[["NGSRELATE_REL"]], width = 1500, height = 500)
ggplot(nr) +
  facet_wrap(vars(data), nrow = 1, strip.position = "left") +
  geom_histogram(aes(x = rab, fill = gsub("_.+", "", data)), bins = 20, show.legend = F) + 
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
  facet_wrap(vars(data), nrow = 1, strip.position = "left") +
  geom_histogram(aes(x = mean_Fa, fill = gsub("_.+", "", data)), bins = 20, show.legend = F) + 
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


# Hap.py
happy_file <- snakemake@input[["HAPPY"]]
happy <- read.csv(happy_file)
happy <- happy %>% 
  pivot_longer(cols = c("METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"),
               names_to = "Metrics", values_to = "Value") 
happy$Metrics <- factor(happy$Metrics, levels = c("METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"))

png(snakemake@output[["HAPPY"]], width = 800, height = 800)
ggplot(happy) +
  facet_grid(rows = vars(Metrics), scales = "free_x") +
  geom_col(aes(y = Value, x = data, fill = gsub("_.+", "", data)), show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1)) + 
  ylab("Value") + 
  xlab("Imputation") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()

# Customized accuracy metrics
var_file <- snakemake@input[["METRICS_VAR"]]
var <- read.delim(var_file)

var <- var %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQS"),
               names_to = "Metrics", values_to = "Value") 
var$Truth_MAF <- sapply(var$Truth_ALT_Freq, function(x) min(c(x, 1-x)))
var$MAF_bin <- cut(var$Truth_MAF, breaks = c(0, 0.05, 0.1, 0.25, 0.5), include.lowest = T)
var$Metrics <- factor(var$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Recall", "Precision", "F1", "IQS"))

png(snakemake@output[["VARIANT"]], width = 1600, height = 900)
ggplot(var) +
  facet_grid(cols = vars(data), rows = vars(Metrics)) +
  geom_boxplot(aes(y = Value, x = MAF_bin, color = gsub("_.+", "", data)), outliers = F, show.legend = F) +
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


sam_file <- snakemake@input[["METRICS_SAM"]]
sam <- read.delim(sam_file)
sam <- sam %>% 
  pivot_longer(cols = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQS"),
               names_to = "Metrics", values_to = "Value") 
sam$Metrics <- factor(sam$Metrics, levels = c("Mismatch_1", "Mismatch_2", "Missing", "Precision", "Recall", "F1", "IQ"))

png(snakemake@output[["SAMPLE"]], width = 1000, height = 800)
ggplot(sam) +
  facet_grid(rows = vars(Metrics), scales = "free_x") +
  geom_boxplot(aes(y = Value, x = data, color = gsub("_.+", "", data)), outliers = F, show.legend = F) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        panel.background = element_rect(color = "black", fill = NULL, linewidth = 1),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1)) + 
  ylab("Value per sample") + 
  xlab("Imputation") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
dev.off()

# MAF
var$Imputed_MAF <- sapply(var$Imputed_ALT_Freq, function(x) min(c(x, 1-x)))

png(snakemake@output[["MAF"]], width = 1500, height = 300)
ggplot(var, aes(y = Imputed_MAF, x = Truth_MAF, color = gsub("_.+", "", data))) +
  facet_wrap(vars(data), nrow = 1, strip.position = "left") +
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
