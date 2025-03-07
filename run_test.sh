## activate snakemake
conda activate snakemake
mamba activate snakemake

## buld docker container
sudo docker build --no-cache -f Dockerfile
sudo docker build -t vibaotram/imputation:1.4 .
sudo docker push vibaotram/imputation:1.4
sudo docker run -it vibaotram/imputation:1.4

## prepare ref vcf
bcftools view -m2 -M2 -v snps -Oz -o test/INPUT/REF_VCF/p1_random_1000_30gen_ref1_filtered.vcf.gz test/INPUT/REF_VCF/p1_random_1000_30gen_ref1.vcf
# bcftools view -Oz -o test/INPUT/REF_VCF/p1_random_1000_30gen_ref1.vcf.gz test/INPUT/REF_VCF/p1_random_1000_30gen_ref1.vcf
bcftools index test/INPUT/REF_VCF/p1_random_1000_30gen_ref1_filtered.vcf.gz
tabix test/INPUT/REF_VCF/p1_random_1000_30gen_ref1_filtered.vcf.gz

## run test
snakemake -n

rm .snakemake/singularity/*
snakemake --sdm apptainer -F
snakemake ACCURACY --sdm apptainer --configfile config.yaml --rerun-triggers mtime -np 
snakemake ACCURACY --sdm apptainer --configfile config.yaml --report report.html

# dag, rulegraph, filegraph, summary
snakemake ACCURACY --configfile config.yaml --rulegraph | dot -Tpdf > test/dag.pdf 


NOTE:

- QUILT doesn't work for small reference panel (e.g. 20 ref with 1400 SNPs after filtering for biallelic sites)
- Customized script: precision = recall in most case where no missingness, except for STITCH that has sporadic missing genotypes