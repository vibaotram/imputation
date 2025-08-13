# Performing imputation and population structure for low-coverage sequencing data, and comparisons with ground truth

The pipeline consists of 3 main parts: imputation by multiple choices of tool, population genetic analysis on imputed genotypes, and evaluation of imputation performance if ground truth is present. The entire pipeline runs on readily available docker images.

## Table of Contents
- [Input data](#input-data)
- [Imputation](#imputation)
- [Population genetic analysis](#population-genetic-analysis)
- [Imputation performance](#imputation-performance)
- [Software requirements](#software-requirement)
- [Run the pipeline](#run-the-pipeline)


## Input data

- TARGET_BAM:       a list of BAM files for imputation targets
- REFERENCE_SEQ:    a reference genome in fasta format
- REFERENCE_VCFGZ:  reference panel in vcf.gz format, index file is required 
- GROUNDTRUTH:      ground truth genotypes in vcf.gz format, if available


## Imputation

Five imputation tools are available and all can be executed in the same pipeline:
- BEAGLE 5.4
- GLIMPSE 2.0.1
- GeneImp 1.3.0.9000
- QUILT 2.0.1
- STITCH 1.7.1

Prior to imputation, reference VCF file will be filtered out for fixed alleles. If BEAGLE or GeneImp is executed, SNP calling and filtering will be performed by bcftools mpileup and bcftools filter, then the output VCF file will be used for imputation.

## Population genetic analysis

Analyses will be performed on the imputed genotypes:
- pairwise relatedness and inbreeding coefficients, by [ANGSD](https://www.popgen.dk/angsd/index.php/ANGSD)/[NgsRelate](https://github.com/ANGSD/NgsRelate) and [KGD](https://github.com/AgResearch/KGD)
- Demographic inference by vcf2sfs and stairway plot v2
- Fst between populations, if more than 1 population are specified

## Evaluation of imputation performance

If a ground truth is present, imputation accuracy will be measured by comparing the imputed genotypes and the truth genotypes. Different metrics, true positives (TP), false positive (FP), false negative (FN), recall, precisiom, F1 score, concordance (Po), chance agreement (Pc), imputation quality score (IQS) will be calculated the basis of variant and sample. 

- TP: number of matching alleles
- FP: number of mismatched alleles and missing alleles
- FN: number of missing alleles
- recall = TP / (TP + FN)
- precision = TP / (TP + FP)
- F1 score = recall * precision / (recall + precision)
- Po, Pc, and IQS are described [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009697)

Imputed genotypes are also compared with the ground truth by [hap.py](https://github.com/Illumina/hap.py) for haplotype comparisons.

In addition, population genetic analyses will also be performed on the ground truth for assessing the effects of imputation on downstream analysis.

## Software requirement
- Snakemake >=8.20.0
- Singularity >= 3.0.0

## Run the pipeline
1. Prepare [input data](#input-data)


Some scripts to test the pipeline are written in this [bash file](./run_test.sh)