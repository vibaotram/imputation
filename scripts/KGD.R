workdir <- dirname(snakemake@output[[1]])
# workdir <- "/home/tvi069/imputation/test/OUTPUT/KGD/IMPUTED/GENEIMP"
setwd(workdir)

genofile <- snakemake@input[[1]]
# genofile <- "/home/tvi069/imputation/test/OUTPUT/GENEIMP/p4_GENEIMP.vcf.gz"
gform <- "VCF"
cachedir <- workdir
source(snakemake@params[[1]])
# source("/home/tvi069/imputation/scripts/KGD/GBS-Chip-Gmatrix.R")
Gfull <- calcG()

# png(snakemake@output[[1]])
# hist(as.vector(Gfull$G5), main = "Relatedness (G5)", xlab = "Relatedness")
# dev.off()

writeLines(as.character(Gfull$G5), snakemake@output[[1]])
writeLines(as.character(Gfull$G4d), snakemake@output[[2]])