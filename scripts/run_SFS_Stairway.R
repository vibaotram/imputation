source("/vcf2sfs/vcf2sfs.r")

vcffile <- snakemake@input[["VCF"]]
popfile <- snakemake@input[["POPFILE"]]
template_bp <- snakemake@params[["BLUEPRINT"]]
outdir <- snakemake@params[['OUTDIR']]
seqlength <- snakemake@params[['SEQ_LEN']]
mu <- snakemake@params[['MU']]
plot_title <- snakemake@param[['PLOT_TITLE']]
output <- snakemake@output[[1]]
stairway_plot_es_path <- "/stairway_plot_v2.1.2/stairway_plot_es"

pops <- read.table(popfile)


gt<-vcf2gt(vcffile, popfile)

blueprint <- readLines(template_bp)
names(blueprint) <- gsub(":.+", "", blueprint)
blueprint["L"] <- paste(c("L:", seqlength), collapse = " ")
blueprint["mu"] <- paste(c("mu:", mu), collapse = " ")
blueprint["project_dir"] <- "project_dir: two-epoch_fold"

for (pop in unique(pops[,2])) {
  sfsraw <- gt2sfs.raw(gt, pop)
  folded_sfs <- fold.sfs(sfsraw)
  nind<- length(pops$V2[pops$V2 == pop])
  
  blueprint["SFS"] <- paste(c("SFS:", folded_sfs[1:nind]), collapse = " ")
  blueprint["popid"] <- paste(c("popid:", pop), collapse = " ")
  blueprint["nseq"] <- paste(c("nseq:", nind*2), collapse = " ")
  blueprint["plot_title"] <- paste(plot_title, "- subpop", pop)
  pop_dir <- paste0(outdir, "/", pop)
  
  if (dir.exists(pop_dir)) unlink(pop_dir, recursive = T)
  dir.create(pop_dir, recursive = T)
  file.copy(stairway_plot_es_path, pop_dir, recursive = T)
  setwd(pop_dir)
  system("pwd")
  updated_bp <- paste0(pop_dir, "/", pop, ".blueprint")
  
  writeLines(blueprint, updated_bp)
  
  stairway_build_cmd <- paste("java -cp stairway_plot_es Stairbuilder", updated_bp)
  system(stairway_build_cmd)
  
  run_staiway_cmd <- paste0("bash ", updated_bp, ".sh")
  system(run_staiway_cmd)
  
  unlink("stairway_plot_es", recursive = T)
}

final_file <- paste0(outdir, "/", unique(pops[,2]), "/", blueprint["project_dir"], "/", blueprint["project_dir"], ".final.summary[.png|.pdf]")
writeLines(final_file, output)
