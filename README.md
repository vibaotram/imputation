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

## Output
### Imputation

Five imputation tools are available and all can be executed in the same pipeline:
- BEAGLE 5.4
- GLIMPSE 2.0.1
- GeneImp 1.3.0.9000
- QUILT 2.0.1
- STITCH 1.7.1

Prior to imputation, reference VCF file will be filtered out for fixed alleles. If BEAGLE or GeneImp is executed, SNP calling and filtering will be performed by bcftools mpileup and bcftools filter, then the output VCF file will be used for imputation.

**Imputed VCF files are created in <OUTPUT_DIR>/<TOOL>/<PARAMETER>/<PREFIX>_<TOOL>.vcf.gz.**

### Population genetic analysis

Analyses will be performed on the imputed genotypes:
- pairwise relatedness and inbreeding coefficients, by [ANGSD](https://www.popgen.dk/angsd/index.php/ANGSD)/[NgsRelate](https://github.com/ANGSD/NgsRelate) and [KGD](https://github.com/AgResearch/KGD)
- Demographic inference by vcf2sfs and stairway plot v2
- Fst between populations, if more than 1 population are specified

**Results are stored in <OUTPUT_DIR>/POPGEN/<ANALYSIS> folders.**

### Evaluation of imputation performance

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

**Output data include:**


## Software requirement
- Snakemake >=8.20.0
- Singularity >= 3.0.0

## Run the pipeline
1. Prepare [input data](#input-data)

2. Configure the pipeline
- Prepare a config file with this [template](config/config.yaml). 
- For each imputation tool - BEAGLE, STITCH, GLIMPSE2 (GLIMPSE2_phase), QUILT2, and GENEIMP, parameters can be provided in a csv file with different parameters in columns. Column names must match exactly parameter name, without the leading "--/-". Different sets of values for the parameters can be provided in different rows. If the csv file is not provided in the config file, default parameters will be used. Details of parameters for all the imputation tools can be found [here](config/params)
- Number of CPUs (THREADS) for each tool in the config file must be specified.

3. Run Snakemake:
The pipeline can run upon 3 main steps (IMPUTATION, POPGEN, and ACCURACY evaluation) and it needs to be specified in the Snakemake command as below.

```
## run imputation only (no downstream analysis)
snakemake IMPUTATION --sdm apptainer --configfile config/config.yaml --cores 12

## run imputation and downstream population genetic analysis
snakemake POPGEN --sdm apptainer --configfile config/config.yaml --cores 12

## run imputation, downstream analysis, and accuracy evaluation if a groundtruth is provided
snakemake ACCURACY --sdm apptainer --configfile config/config.yaml --cores 12
```
A test data set of 10 LCS samples, with 10 samples in the reference panel, and a groundtruth set of SNPs, is provided in the `test` folder. Scripts to run the test can be found in this [bash file](./run_test.sh).

## Note
Snakemake sometimes does not automatically detect the total number of cores for the pipeline, and users need to specified "--cores <N>".