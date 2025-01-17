#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
library(GeneImp)

option_list <- list(
  make_option(
    "--target",
    type = "character",
    help = "Path to the VCF file containing the genotype likelihoods for the target samples. Genotype likelihoods are read from the PL field. The file must contain reads from a single chromosome."
  ), 
  make_option(
    "--ref",
    type = "character",
    help = "Path to the VCF file containing phased genotypes for the reference panel. Reference haplotypes are read from the GT field which is always assumed to contain phased genotypes. The file must contain variants from a single chromosome."
  ), 
  make_option(
    "--klthresh",
    type = "numeric",
    default = 20,
    help = "Size of the imputation window. Larger values will result in larger windows. The klthresh option determines the size of the imputation window. In target data with average sequencing depth 0.5x setting the klthresh to 15, 20, 25, 30, 35 resulted in windows with median length 58, 79, 100, 121, 142kb respectively. In data with deeper sequencing smaller values for the klthresh should be used to get windows with similar median lengths. Only the central part of each window is imputed based on the flanksize option."
  ), 
  make_option(
    "--flanksize",
    type = "numeric",
    help = "	
Size of the flanking segments left and right of the imputation window. The flanksize option determines the size of the overlap between windows). Only the central part of each window is imputed. A value of 0.5 means that approximately 25% of the reads come from the left flanking segment, 50% of the reads from the window to be imputed, and 25% of the reads from the right flanking segment.",
    default = 0.5
  ), 
  make_option(
    "--filtermethod",
    type = "character",
    default = "pairrand",
    help = "Method for selecting reference haplotypes. Choose one of these: 'pairrand', 'paircond', 'marginal', 'random', or 'none'"
  ), 
  make_option(
    "--numfilterhaps",
    type = "integer",
    default = 20,
    help = "Number of reference haplotypes to select in the filtering step."
  ), 
  make_option(
    "--maxjobs",
    type = "integer",
    default = 1,
    help = "Maximum number of parallel processes to run with doMC."
  ), 
  make_option(
    "--writedir",
    type = "character",
    help = "Directory where the imputed data will be saved. Specify if different from current directory.",
    default = ""
  ), 
  make_option(
    "--tempdir",
    type = "character",
    help = "Directory for storing large files used during the imputation, such as the file.backed bigmatrix containing the reference haplotypes. Specify if different from current directory.",
    default = ""
  ), 
  make_option(
    "--verbose",
    type = "integer",
    default = 1,
    help = "Print various messages as the code executes."
  ), 
  make_option(
    "--diagnostics",
    type = "logical",
    default = FALSE,
    help = "Set to TRUE to write imputation diagnostics, such as window positions, filtered haplotypes and haplotype probabilities. Files are written in the temp.dir if specified. CAREFUL: these files are written for every sample."
  )
)


opt <- parse_args(OptionParser(option_list = option_list))

imputevcf(vcfname = opt$target, 
          ref.vcfname = opt$ref, 
          klthresh = opt$klthresh, 
          flanksize = opt$flanksize, 
          filtermethod = opt$filtermethod,
          numfilterhaps = opt$numfilterhaps, 
          maxjobs = opt$maxjobs, 
          write.dir = opt$writedir, 
          temp.dir = opt$tempdir, 
          verbose = opt$verbose,
          diagnostics = opt$diagnostics)