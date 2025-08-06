## activate snakemake
conda activate snakemake
mamba activate snakemake

## run imputation only
snakemake IMPUTATION --sdm apptainer --configfile config/config.yaml

## run imputation and assessment if a ground truth is present
snakemake ACCURACY --sdm apptainer --configfile config/config.yaml

## generate a report file containing all the output figures
snakemake ACCURACY --sdm apptainer --configfile config/config.yaml --report test/OUTPUT/report.html

## create dag, rulegraph, filegraph, or summary
snakemake ACCURACY --configfile configconfig.yaml --dag [rulegraph,filegraph,summary] | dot -Tpdf > test/dag.pdf 
