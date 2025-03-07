library(vcfR)
library(parallel)
library(data.table)

truth <- snakemake@input[["TRUTH"]]
imputed <- snakemake@input[["TEST"]]
ncores <- snakemake@threads

truth_data <- read.vcfR(truth)
imputed_data <- read.vcfR(imputed)

truth_info <- getFIX(truth_data)
truth_snp <- paste(truth_info[,1], truth_info[,2], sep = "_")
truth_gt <- extract.gt(truth_data)
rownames(truth_gt) <- truth_snp

imputed_info <- getFIX(imputed_data)
imputed_snp <- paste(imputed_info[,1], imputed_info[,2], sep = "_")
imputed_gt <- extract.gt(imputed_data)
rownames(imputed_gt) <- imputed_snp



gt_to_num <- function(x) {
  x[x %in% c("0/0", "0|0")] <- 0
  x[x %in% c("0/1", "0|1", "1/0", "1|0")] <- 1
  x[x %in% c("1/1", "1|1")] <- 2
  x <- apply(x, c(1,2), as.numeric)
  return(x)
}




accuracy_metrics <- function(actual, predicted, ncores) {
  # actual <- truth_gt
  # predicted <- imputed_gt
  
  cat("Converting gt (0/0, 0/1, 1/1, ./.) into number (0,1,2, NA)\n")
  # convert gt (0/0, 0/1, 1/1, ./.) into number (0,1,2, NA)
  actual_num <- gt_to_num(actual)
  predicted_num <- gt_to_num(predicted)
  common_ind <- colnames(actual_num)[colnames(actual_num) %in% colnames(predicted_num)]
  actual_num <- actual_num[, colnames(actual_num) %in% common_ind]
  predicted_num <- predicted_num[, colnames(predicted_num) %in% common_ind]

  # fill missing sites in predicted matrix by NA 
  na_mat <- actual_num
  na_mat[!is.na(na_mat)] <- NA
  cat("predicted_num:", dim(predicted_num), "\n")
  cat("na_mat:", dim(actual_num), "\n")
  predicted_num <- rbind(predicted_num, na_mat[!rownames(na_mat) %in% rownames(predicted_num),])
  predicted_num <- predicted_num[rownames(actual_num), colnames(actual_num)]
  
  cat("Counting TP, FP, FN, and mismatches per variant and per sample. This might take a while\n")
  counts <- mclapply(rownames(actual_num), function(i) {
    mapply(function(x, y, TP = NA, FP = NA, FN = NA, mismatch_1 = NA, mismatch_2 = NA, missing = NA, abs_match = 0) {
      if (!is.na(x) && !is.na(y)) { ## if no missing gt in both sets, counts the number of matched and mismatched
        mismatched <- abs(x - y)
        matched <- 2 - mismatched
        TP <- matched
        ## mismatched is counted as both FP and FN
        FP <- mismatched
        FN <- mismatched
        if (mismatched == 1) {
          mismatch_1 <- 1
          mismatch_2 <- 0
        } else if (mismatched == 2) {
          mismatch_1 <- 0
          mismatch_2 <- 1
        } else {
          mismatch_1 <- 0
          mismatch_2 <- 0
          abs_match <- 1
        }
        missing <- 0
      } else if (!is.na(x) && is.na(y)) { ## if predicted gt is missing, count 2 for FN
        TP <- 0
        FP <- 0
        FN <- 2
        mismatch_1 <- 0
        mismatch_2 <- 0
        missing <- 2
      } ## if actual gt is missing, do nothing
      return(c(TP, FP, FN, mismatch_1, mismatch_2, missing, abs_match))
    }, actual_num[i,], predicted_num[i,])
  }, mc.cores = ncores)
  
  cat("Calculating Recall, Precision, F1 score, IQS\n")
  actual_alt_freq <- rowMeans(actual_num, na.rm = T)/2
  predicted_alt_freq <- rowMeans(predicted_num, na.rm = T)/2
  
  
  TP_mat <- t(sapply(counts, function(c) as.numeric(c[1,])))
  colnames(TP_mat) <- colnames(actual_num)
  TP_per_variant <- rowSums(TP_mat, na.rm = T)
  TP_per_sample <- colSums(TP_mat, na.rm = T)
  
  
  FP_mat <- t(sapply(counts, function(c) as.numeric(c[2,])))
  colnames(FP_mat) <- colnames(actual_num)
  FP_per_variant <- rowSums(FP_mat, na.rm = T)
  FP_per_sample <- colSums(FP_mat, na.rm = T)
  
  
  FN_mat <- t(sapply(counts, function(c) as.numeric(c[3,])))
  colnames(FN_mat) <- colnames(actual_num)
  FN_per_variant <- rowSums(FN_mat, na.rm = T)
  FN_per_sample <- colSums(FN_mat, na.rm = T)
  
  
  mm1_mat <- t(sapply(counts, function(c) as.numeric(c[4,])))
  colnames(mm1_mat) <- colnames(actual_num)
  mm1_per_variant <- rowMeans(mm1_mat, na.rm = T)/2
  mm1_per_sample <- colMeans(mm1_mat, na.rm = T)/2
  
  
  mm2_mat <- t(sapply(counts, function(c) as.numeric(c[5,])))
  colnames(mm2_mat) <- colnames(actual_num)
  mm2_per_variant <- rowMeans(mm2_mat, na.rm = T)/2
  mm2_per_sample <- colMeans(mm2_mat, na.rm = T)/2
  
  miss_mat <- t(sapply(counts, function(c) as.numeric(c[6,])))
  colnames(miss_mat) <- colnames(actual_num)
  miss_per_variant <- rowMeans(miss_mat, na.rm = T)/2
  miss_per_sample <- colMeans(miss_mat, na.rm = T)/2
  
  absmatch_mat <- t(sapply(counts, function(c) as.numeric(c[7,])))
  colnames(absmatch_mat) <- colnames(actual_num)
  po_per_variant <- rowMeans(absmatch_mat, na.rm = T)
  po_per_sample <- colMeans(absmatch_mat, na.rm = T)
  
  pc_per_variant <- sapply(1:nrow(actual_num), function(r) {
    freq_actual <- (table(c(actual_num[r,], 0,1,2))-1)/length(actual_num[r,])
    freq_predicted <- (table(c(predicted_num[r,], 0,1,2))-1)/length(actual_num[r,])
    pc <- sum(freq_actual*freq_predicted)
    pc
  })
  iqs_per_variant <- (po_per_variant-pc_per_variant)/(1-pc_per_variant)
  
  pc_per_sample <- sapply(1:ncol(actual_num), function(c) {
    freq_actual <- (table(c(actual_num[,c], 0,1,2))-1)/length(actual_num[,c])
    freq_predicted <- (table(c(predicted_num[,c], 0,1,2))-1)/length(actual_num[,c])
    pc <- sum(freq_actual*freq_predicted)
    pc
  })
  iqs_per_sample <- (po_per_sample-pc_per_sample)/(1-pc_per_sample)
  
  var <- data.frame(ID = rownames(actual_num),
                    Truth_ALT_Freq = actual_alt_freq,
                    Imputed_ALT_Freq = predicted_alt_freq,
                    Mismatch_1 = mm1_per_variant,
                    Mismatch_2 = mm2_per_variant,
                    Missing = miss_per_variant,
                    Po = po_per_variant,
                    Pc = pc_per_variant,
                    IQS = iqs_per_variant,
                    TP = TP_per_variant,
                    FP = FP_per_variant,
                    FN = FN_per_variant)
  var$Recall <- var$TP/(var$TP + var$FN)
  var$Precision <- var$TP/(var$TP + var$FP)
  var$F1 <- 2 * var$Recall * var$Precision/ (var$Recall + var$Precision)
  
  
  sam <- data.frame(ID = colnames(actual_num),
                    Mismatch_1 = mm1_per_sample,
                    Mismatch_2 = mm2_per_sample,
                    Missing = miss_per_sample,
                    Po = po_per_sample,
                    Pc = pc_per_sample,
                    IQS = iqs_per_sample,
                    TP = TP_per_sample,
                    FP = FP_per_sample,
                    FN = FN_per_sample)
  sam$Recall <- sam$TP/(sam$TP + sam$FN)
  sam$Precision <- sam$TP/(sam$TP + sam$FP)
  sam$F1 <- 2 * sam$Recall * sam$Precision / (sam$Recall + sam$Precision)
  
  cat("Finished!")
  return(list("variant" = var, "sample" = sam))
}

metrics <- accuracy_metrics(truth_gt, imputed_gt, ncores)

fwrite(metrics$variant, snakemake@output[[1]], sep = "\t")
fwrite(metrics$sample, snakemake@output[[2]], sep = "\t")

# imputed <- "~/imputation/test/OUTPUT/ACCURACY/BEAGLE/p4_GT.vcf.gz"
# imputed <- "~/imputation/test/OUTPUT/ACCURACY/STITCH/p4_GT.vcf.gz"
# truth <- "~/imputation/test/OUTPUT/VCF/p4_GROUNDTRUTH.vcf.gz"
# t <- read.table("~/imputation/test/OUTPUT/ACCURACY/BEAGLE/p4_BEAGLE_per_variant_results.txt", header = T)
# t <- read.table("~/imputation/test/OUTPUT/ACCURACY/STITCH/p4_STITCH_per_variant_results.txt", header = T)
# t$position <- as.numeric(gsub(".+:", "", t$position))

# metrics$variant$position <- rownames(metrics$variant)
# metrics$variant$position <- gsub(".+\\_", "", metrics$variant$position)

# d <- merge(t, metrics$variant, by = "position")

# ggplot(d, aes(x = IQS.x, y = IQS.y)) + geom_point() + theme_minimal() +
#   xlab("IQS by github tool") +
#   ylab("IQS by custom script")

# metrics_variant %>% 
#   filter(!is.nan(Imputed_ALT_Freq)) %>% 
#   pivot_longer(cols = Recall:F1, names_to = "Metrics", values_to = "Value") %>%
#   ggplot() +
#   # facet_grid(rows = vars(Metrics)) + 
#   geom_boxplot(aes(x = Metrics, y = Value), outliers = F,
#                position = position_dodge2(width = 0.5, padding = 0))
