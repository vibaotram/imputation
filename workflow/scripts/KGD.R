workdir <- dirname(snakemake@output[[1]])

setwd(workdir)

genofile <- snakemake@input[[1]]

gform <- "VCF"
cachedir <- workdir
source(snakemake@params[[1]])

Gfull <- calcG()


writeLines(as.character(Gfull$G5), snakemake@output[[1]])
writeLines(as.character(Gfull$G4d), snakemake@output[[2]])