workdir <- dirname(snakemake@output[[1]])
setwd(workdir)

genofile <- snakemake@input[[1]]
gform <- "VCF"
cachedir <- workdir
source(snakemake@params[[1]])
Gfull <- calcG()

# png(snakemake@output[[1]])
# hist(as.vector(Gfull$G5), main = "Relatedness (G5)", xlab = "Relatedness")
# dev.off()

writeLines(as.character(Gfull$G5), snakemake@output[[1]])