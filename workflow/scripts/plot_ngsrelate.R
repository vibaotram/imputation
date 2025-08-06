library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ngsrelate_file <- args[1]

ngsrelate <- read.delim(ngsrelate_file)

png(args[2])
hist(ngsrelate$rab, main = "Pairwise Relatedness (rab)", xlab = "rab")
dev.off()


ngsrelate_fa <- ngsrelate %>% 
  group_by(a) %>% 
  summarise(mean_Fa = mean(Fa))


png(args[3])
hist(ngsrelate_fa$mean_Fa, main = "Individual inbreeding coefficient (Fa)", xlab = "Fa")
dev.off()
