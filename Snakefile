import pandas as pd
container: "docker://vibaotram/imputation:1.4" #docker pull staphb/bcftools, continuumio/miniconda3:24.7.1-0

#configfile: "config.yaml"


IMPUTATION = []
if config['IMPUTATION']['BEAGLE']:
    IMPUTATION.append("BEAGLE")
if config['IMPUTATION']['GLIMPSE']:
    IMPUTATION.append("GLIMPSE")
if config['IMPUTATION']['STITCH']:
    IMPUTATION.append("STITCH")
if config['IMPUTATION']['GENEIMP']:
    IMPUTATION.append("GENEIMP")
if config['IMPUTATION']['QUILT']:
    IMPUTATION.append("QUILT")


with open(config["TARGET_BAM"]) as BAMlist:
    BAMfiles = []
    for line in BAMlist:
        BAMfiles.append(line)


CHROM = []
REFfile = open(config['REFERENCE_SEQ'])
for line in REFfile:
    if ">" in line:
        CHROM.append(line[1:-1])
REFfile.close()

rule IMPUTATION:
    input:
        expand(os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz"), tool = IMPUTATION)


#-------------PREPARE--------------#
#-------------TARGET---------------#

rule NGSRELATE_BAM:
    input: config['TARGET_BAM']
    output:
        RESULTS = os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "BAM", config['PREFIX_NAME'] + ".ngsrelate")
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "BAM"),
        OUTNAME = config['PREFIX_NAME'],
        NIND = len(BAMfiles),
        MIN_MAF = config['NGSRELATE']['MIN_MAF']
    threads: config['NGSRELATE']['THREADS']
    # conda: "envs/angsd.yaml"
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "NGSRELATE.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        mkdir -p {params.OUTDIR}

        # genotype likelihood and allele freq
        angsd -b {input} -gl 2 \
        -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf {params.MIN_MAF} -doGlf 3 \
        -out {params.OUTDIR}/{params.OUTNAME} \
        -P {threads}

        zcat {params.OUTDIR}/{params.OUTNAME}.mafs.gz | cut -f5 | sed 1d > {params.OUTDIR}/{params.OUTNAME}.freq
        
        id_list=$(mktemp -p {params.OUTDIR})
        while read FILE;
        do
            basename $FILE | sed 's/.bam//' >> $id_list
        done < {input}

        ### run NgsRelate
        ngsRelate -g {params.OUTDIR}/{params.OUTNAME}.glf.gz \
        -n {params.NIND} \
        -z $id_list \
        -f {params.OUTDIR}/{params.OUTNAME}.freq  \
        -O {output.RESULTS} \
        -p {threads}
        
        rm $id_list
        '''

rule BAM2VCF:
    input:
        BAMlist = config['TARGET_BAM'],
        REF_SEQ = config['REFERENCE_SEQ']
    output: os.path.join(config['OUTPUT_DIR'], "VCF", config['PREFIX_NAME'] + ".vcf.gz")
    threads: config['BAM2VCF']['THREADS']
    # conda: "envs/bcftools.yaml"
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "BAM2VCF.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        bcftools mpileup -a 'DP,DP4,AD,ADF,ADR,SP' --threads {threads} -Ou \
        -f {input.REF_SEQ} \
        -b {input.BAMlist} | \
        bcftools call -mv -G - | \
        bcftools view -m2 -M2 -v snps | \
        bcftools sort --temp-dir $(dirname {output})/tmp -Oz -o {output}
        '''


rule FILTERVCF_REF:
    input:
        REF = config['REFERENCE_VCFGZ'],
        TRUTH = config['GROUNDTRUTH'] if config['GROUNDTRUTH'] else ()
    output:
        REF = os.path.join(config['OUTPUT_DIR'], "VCF", config['PREFIX_NAME'] + "_REF.vcf.gz"),
        TRUTH = os.path.join(config['OUTPUT_DIR'], "VCF", config['PREFIX_NAME'] + "_GROUNDTRUTH.vcf.gz") if config['GROUNDTRUTH'] else ()
    # params:
    #     MIN_QUAL = config['FILTERVCF']['MIN_QUAL']
    threads: config['FILTERVCF_TARGET']['THREADS']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "FILTERVCF_REF.log")
    shell:
        '''
        # filter REF for fixed alleles
        bcftools view -c 1:minor --threads {threads} -Oz -o {output.REF} {input.REF}
        tabix -f {output.REF}

        # filter GROUNDTRUTH for same sites as in REF
        if [[ -n {input.TRUTH} ]]; then
            bcftools view -R {output.REF} --threads {threads} -Oz -o {output.TRUTH} {input.TRUTH}
            tabix -f {output.TRUTH}
        fi
        '''



rule FILTERVCF_TARGET:
    input:
        TARGET = rules.BAM2VCF.output,
        REF = rules.FILTERVCF_REF.output.REF
    output: os.path.join(config['OUTPUT_DIR'], "VCF", config['PREFIX_NAME'] + "_filtered.vcf.gz")
    params:
        MIN_QUAL = config['FILTERVCF_TARGET']['MIN_QUAL']
    threads: config['FILTERVCF_TARGET']['THREADS']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "FILTERVCF_TARGET.log")
    shell:
        '''
        # filter TARGET
        bcftools index -f {input.TARGET}
        bcftools filter -R {input.REF} -i 'QUAL >= {params.MIN_QUAL}' --threads {threads} -Oz -o {output} {input.TARGET}
        bcftools index -f {output}
        tabix -f {output}
        '''


rule KGD_LCS:
    input: rules.FILTERVCF_TARGET.output
    output:
        RELATEDNESS = os.path.join(config['OUTPUT_DIR'], "KGD", "LCS", config['PREFIX_NAME'] + "_pairwise_relatedness.txt"),
        INBREEDING = os.path.join(config['OUTPUT_DIR'], "KGD", "LCS", config['PREFIX_NAME'] + "_inbreeding_coefficient.txt")
    params:
        KGD_GMATRIX = os.path.abspath("scripts/KGD/GBS-Chip-Gmatrix.R")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "KGD.log")
    script:
        "scripts/KGD.R"

rule KGD_TRUTH:
    input: rules.FILTERVCF_REF.output.TRUTH
    output:
        RELATEDNESS = os.path.join(config['OUTPUT_DIR'], "KGD", "TRUTH", config['PREFIX_NAME'] + "_pairwise_relatedness.txt"),
        INBREEDING = os.path.join(config['OUTPUT_DIR'], "KGD", "TRUTH", config['PREFIX_NAME'] + "_inbreeding_coefficient.txt")
    params:
        KGD_GMATRIX = os.path.abspath("scripts/KGD/GBS-Chip-Gmatrix.R")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "KGD.log")
    script:
        "scripts/KGD.R"


rule NGSRELATE_TRUTH:
    input: rules.FILTERVCF_REF.output.TRUTH
    output:
        RESULTS = os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "TRUTH", config['PREFIX_NAME'] + ".ngsrelate")
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "NGSRELATE"),
        OUTNAME = config['PREFIX_NAME'],
        NIND = len(BAMfiles),
        MIN_MAF = config['NGSRELATE']['MIN_MAF']
    threads: config['NGSRELATE']['THREADS']
    # conda: "envs/angsd.yaml"
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "NGSRELATE.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1

        id_list=$(mktemp -p {params.OUTDIR})
        bcftools query -l {input} > $id_list
        ngsRelate -h {input} -T GT -z $id_list -O {output.RESULTS} -p {threads}
        rm $id_list
        '''


#--------------PREPARE-------------#
#-------------REFERENCE------------#

rule REF_BEAGLE_PHASING:
    input: rules.FILTERVCF_REF.output.REF # config['REFERENCE_VCFGZ']
    output: os.path.join(config['OUTPUT_DIR'], "PHASE", config['PREFIX_NAME'] + "_phased.vcf.gz")
    threads: config['REF_BEAGLE_PHASING']['THREADS']
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "PHASE"),
        TMP_OUTNAME = config['PREFIX_NAME'] + "_TMP",
        PARAM = config['REF_BEAGLE_PHASING']['PARAM']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "REF_BEAGLE_PHASING.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        java -jar /beagle-5.4.jar \
        gt={input} \
        out={params.OUTDIR}/{params.TMP_OUTNAME} \
        {params.PARAM}
        tabix -f {params.OUTDIR}/{params.TMP_OUTNAME}.vcf.gz

        bcftools +fill-tags {params.OUTDIR}/{params.TMP_OUTNAME}.vcf.gz --threads {threads} -Oz -o {output} -- -t AN,AC
        tabix -f {output}
        rm {params.OUTDIR}/{params.TMP_OUTNAME}.vcf.gz
        '''

#----------------RUN---------------#
#--------------BEAGLE--------------#

rule BEAGLE:
    input:
        TARGET = rules.FILTERVCF_TARGET.output,
        REF = rules.REF_BEAGLE_PHASING.output
    output: os.path.join(config['OUTPUT_DIR'], "BEAGLE", config['PREFIX_NAME'] + "_BEAGLE.vcf.gz")
    threads: config['BEAGLE']['THREADS']
    params:
        OUT = os.path.join(config['OUTPUT_DIR'], "BEAGLE", config['PREFIX_NAME'] + "_BEAGLE"),
        #NE = config['BEAGLE']['NE'],
        #PHASE_STATES = config['BEAGLE']['PHASE_STATES'],
        PARAM = config['BEAGLE']['PARAM']
    benchmark: os.path.join(config['OUTPUT_DIR'], "BENCHMARK", "BEAGLE.txt")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "BEAGLE.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        java -jar /beagle-5.4.jar \
        ref={input.REF} \
        gt={input.TARGET} \
        out={params.OUT} \
        nthreads={threads} \
        {params.PARAM}

        tabix -f {output}
        '''


#----------------RUN---------------#
#--------------GLIMPSE-------------#

rule GLIMPSE_CHUNK:
    input: rules.REF_BEAGLE_PHASING.output
    output:
        os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "CHUNK", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}_CHUNKS.txt")
    params:
        PARAM = config['GLIMPSE']['GLIMPSE_CHUNK_PARAM'],
        LOG = os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "CHUNK", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}_CHUNKS.log")
    threads: config['GLIMPSE']['THREADS']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GLIMPSE_CHUNK_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        rm -f {output}
        # bcftools view -Oz -o {input} {input}
        # bcftools index {input}
        GLIMPSE2_chunk_static --input {input} --region {wildcards.CHROM} --output {output} --sequential 1 --threads {threads} {params.PARAM} > {params.LOG}

        if [ ! -s {output} ]; then
          REG=$(grep "\\-buffer:" {params.LOG} | cut -d "[" -f3 | cut -d "]" -f1)
          CHR=$(echo $REG | cut -d ":" -f1)
          echo -e "0\t${{CHR}}\t${{REG}}\t${{REG}}" > {output}
        fi
        '''

rule GLIMPSE_SPLIT_REFERENCE:
    input:
        CHUNK = rules.GLIMPSE_CHUNK.output,
        REF_VCF = rules.REF_BEAGLE_PHASING.output
    output: os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "CHUNK", "{CHROM}", config['PREFIX_NAME'] + ".BIN.DONE")
    params:
        BIN = os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "CHUNK", "{CHROM}", config['PREFIX_NAME']),
        PARAM = config['GLIMPSE']['GLIMPSE_SPLIT_REFERENCE_PARAM']
    threads: config['GLIMPSE']['THREADS']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GLIMPSE_SPLIT_REFERENCE_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        rm -f $(dirname {output})/*.bin
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            GLIMPSE2_split_reference_static --reference {input.REF_VCF} --input-region $IRG --output-region $ORG --output {params.BIN} --threads {threads} {params.PARAM}
        done < {input.CHUNK}
        touch {output}
        '''


rule GLIMPSE_PHASE:
    input:
        BAMlist = config['TARGET_BAM'],
        CHUNK = rules.GLIMPSE_CHUNK.output,
        BIN = rules.GLIMPSE_SPLIT_REFERENCE.output
    output: os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "IMPUTE", "{CHROM}", config['PREFIX_NAME'] + ".IMPUTEDLIST")
    params:
        BIN = os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "CHUNK", "{CHROM}", config['PREFIX_NAME']),
        OUTNAME = os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "IMPUTE", "{CHROM}", config['PREFIX_NAME']),
        PARAM = config['GLIMPSE']['GLIMPSE_PHASE_PARAM']
    threads: config['GLIMPSE']['THREADS']
    benchmark: os.path.join(config['OUTPUT_DIR'], "BENCHMARK", "GLIMPSE_PHASE_CHROM{CHROM}.txt")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GLIMPSE_PHASE_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        rm -f {params.OUTNAME}*.bcf
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            CHR=$(echo $LINE | cut -d" " -f2)
            REGS=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f1)
            REGE=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f2)
            GLIMPSE2_phase_static \
            --bam-list {input.BAMlist} \
            --reference {params.BIN}_"$CHR"_"$REGS"_"$REGE".bin \
            --output {params.OUTNAME}_"$REGS"_"$REGE".bcf \
            --threads {threads} {params.PARAM}
        done < {input.CHUNK}
        ls -1v {params.OUTNAME}*.bcf > {output}
        '''


rule GLIMPSE_LIGATE:
    input: rules.GLIMPSE_PHASE.output
    output: os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "LIGATE", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}" + "_GLIMPSE.bcf")
    params:
        PARAM = config['GLIMPSE']['GLIMPSE_LIGATE_PARAM']
    threads: config['GLIMPSE']['THREADS']
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GLIMPSE_LIGATE_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        GLIMPSE2_ligate_static --input {input} --output {output} --threads {threads} {params.PARAM}
        '''

rule GLIMPSE_CONCAT:
    input: expand(os.path.join(config['OUTPUT_DIR'], "GLIMPSE", "LIGATE", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}" + "_GLIMPSE.bcf"), CHROM = CHROM)
    output: os.path.join(config['OUTPUT_DIR'], "GLIMPSE", config['PREFIX_NAME'] + "_GLIMPSE.vcf.gz")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GLIMPSE_CONCAT.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        chr=({CHROM})
        n_chr=${{#chr[@]}}
        if [[ ${{n_chr}} -eq 1 ]]
        then
            bcftools view -Oz -o {output} {input}
        else
            bcftools concat -Oz -o {output} {input}
        fi

        tabix -f {output}
        '''




#----------------RUN---------------#
#--------------STITCH--------------#

rule STITCH:
    input:
        BAMlist = config['TARGET_BAM'],
        REF_VCF = rules.REF_BEAGLE_PHASING.output
    output:
        IMPUTED_VCF = os.path.join(config['OUTPUT_DIR'], "STITCH", "CHROM", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}" + "_STITCH.vcf.gz")
    params:
        TEMPDIR = os.path.join(config['OUTPUT_DIR'], "STITCH", "TEMP", "{CHROM}"),
        POSFILE = os.path.join(config['OUTPUT_DIR'], "STITCH", "TEMP", "{CHROM}", "{CHROM}" + ".POS"),
        OUTNAME = config['PREFIX_NAME'] + "_{CHROM}" + "_STITCH.vcf.gz",
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "STITCH", "CHROM", "{CHROM}"),
        NHAP = config['STITCH']['NHAP'],
        NGEN = config['STITCH']['NGEN'],
        OTHER = config['STITCH']['OTHER']
    threads: config['STITCH']['THREADS']
    benchmark: os.path.join(config['OUTPUT_DIR'], "BENCHMARK", "STITCH_CHROM{CHROM}.txt")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "STITCH_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        # bcftools view -Oz -o {input.REF_VCF} {input.REF_VCF}
        # bcftools index {input.REF_VCF}
        mkdir -p {params.TEMPDIR}
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" -r {wildcards.CHROM} {input.REF_VCF} > {params.POSFILE}

        STITCH --chr={wildcards.CHROM} \
        --bamlist={input.BAMlist} \
        --posfile={params.POSFILE} \
        --outputdir={params.OUTDIR} \
        --output_filename={params.OUTNAME} \
        --tempdir={params.TEMPDIR} \
        --K={params.NHAP} \
        --nGen={params.NGEN} \
        --nCores={threads} \
        {params.OTHER}

        tabix -f {output.IMPUTED_VCF}
        '''


rule STITCH_CONCAT:
    input: expand(rules.STITCH.output.IMPUTED_VCF, CHROM = CHROM)
    output: os.path.join(config['OUTPUT_DIR'], "STITCH", config['PREFIX_NAME'] + "_STITCH.vcf.gz")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "STITCH_CONCAT.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        chr=({CHROM})
        n_chr=${{#chr[@]}}
        if [[ ${{n_chr}} -eq 1 ]]
        then
            cp {input} {output}
        else
            bcftools concat -Oz -o {output} {input}
        fi

        tabix -f {output}
        '''



#----------------RUN---------------#
#--------------GENEIMP-------------#

rule GENEIMP:
    input:
        TARGET = rules.FILTERVCF_TARGET.output,
        REF = rules.REF_BEAGLE_PHASING.output
    output: os.path.join(config['OUTPUT_DIR'], "GENEIMP", config['PREFIX_NAME'] + "_GENEIMP.vcf.gz")
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "GENEIMP"),
        OUTNAME = os.path.basename(rules.FILTERVCF_TARGET.output[0]).split(".")[0],
        KL = config['GENEIMP']['KL'],
        FLANK = config['GENEIMP']['FLANKSIZE'],
        FL = int(float(config['GENEIMP']['FLANKSIZE'])*100),
        HAP = config['GENEIMP']['HAP'],
        FILTER = config['GENEIMP']['FILTER'],
        OTHER = config['GENEIMP']['OTHER']
    threads: config['GENEIMP']['THREADS']
    benchmark: os.path.join(config['OUTPUT_DIR'], "BENCHMARK", "GENEIMP.txt")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "GENEIMP.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        rm -rf {params.OUTDIR}/indivs
        tempdir=$(mktemp -d -p {params.OUTDIR} --suffix {params.OUTNAME})

        Rscript scripts/GeneImp.R \
        --target {input.TARGET} \
        --ref {input.REF} \
        --klthresh {params.KL} \
        --flanksize {params.FLANK} \
        --numfilterhaps {params.HAP} \
        --filtermethod '{params.FILTER}' \
        --maxjobs {threads} \
        --writedir {params.OUTDIR} \
        --tempdir ${{tempdir}} \
        {params.OTHER}

        mv {params.OUTDIR}/all.imputed.{params.OUTNAME}.kl{params.KL}.flank{params.FL}.haps{params.HAP}.{params.FILTER}.vcf.gz {output}
        tabix -f {output}
        '''


#----------------RUN---------------#
#---------------QUILT2-------------#

rule QUILT:
    input:
        BAMlist = config['TARGET_BAM'],
        REF_VCF = rules.REF_BEAGLE_PHASING.output
    output:
        IMPUTED_VCF = os.path.join(config['OUTPUT_DIR'], "QUILT", "CHROM", "{CHROM}", config['PREFIX_NAME'] + "_{CHROM}" + "_QUILT.vcf.gz")
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "QUILT", "CHROM", "{CHROM}"),
        OUTNAME=config['PREFIX_NAME'] + "_{CHROM}" + "_QUILT",
        WINDOW_SIZE = config['QUILT']['WINDOW_SIZE'],
        BUFFER_SIZE = config['QUILT']['BUFFER_SIZE'],
        NGEN = config['QUILT']['NGEN'],
        OTHER = config['QUILT']['OTHER']
    threads: config['QUILT']['THREADS']
    # resources:
        # mem_mb = config['QUILT']['MEM_MB']
    benchmark: os.path.join(config['OUTPUT_DIR'], "BENCHMARK", "QUILT_CHROM{CHROM}.txt")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "QUILT_CHROM{CHROM}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        rm -rf {params.OUTDIR}
        mkdir -p {params.OUTDIR}
        CHR_LEN=$(bcftools view -r {wildcards.CHROM} -Ov {input.REF_VCF} | tail -n 1 | cut -f2)
        SIZE={params.WINDOW_SIZE}
        result=$(echo "scale=2; ${{CHR_LEN}} / ${{SIZE}}" | bc)
        CHUNKS=$(printf "%.0f" "${{result}}")
        echo "Imputing Chunks based on WINDOW_SIZE=$SIZE"
        CHUNK_VCFS=""
        for c in $(seq 1 ${{CHUNKS}})
        do
            start=$(((c - 1) * SIZE + 1))
            if [ ${{c}} -lt ${{CHUNKS}} ]
            then
                end=$((c*SIZE))
            else
                end=${{CHR_LEN}}
            fi
            echo "CHUNK ${{c}}: ${{start}}-${{end}}"

            /QUILT/QUILT2.R \
            --outputdir={params.OUTDIR} \
            --output_filename={params.OUTNAME}_CHUNK${{c}}.vcf.gz \
            --chr={wildcards.CHROM} \
            --regionStart=${{start}} \
            --regionEnd=${{end}} \
            --buffer={params.BUFFER_SIZE} \
            --nGen={params.NGEN} \
            --bamlist={input.BAMlist} \
            --reference_vcf_file={input.REF_VCF} \
            --nCores={threads} \
            {params.OTHER} || \
            echo -e "\e[31mSeems like Quilt2.R has failed either because of small reference panel or memory issues.\nIncreasing memory resource and decreasing threads might fix the memeory issues.\e[0m"


            CHUNK_VCFS="${{CHUNK_VCFS}} {params.OUTDIR}/{params.OUTNAME}_CHUNK${{c}}.vcf.gz"
        done

        bcftools concat --threads {threads} -Oz -o {output.IMPUTED_VCF} ${{CHUNK_VCFS}}
        '''

rule QUILT_CONCAT:
    input: expand(rules.QUILT.output.IMPUTED_VCF, CHROM = CHROM)
    output: os.path.join(config['OUTPUT_DIR'], "QUILT", config['PREFIX_NAME'] + "_QUILT.vcf.gz")
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "QUILT_CONCAT.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        chr=({CHROM})
        n_chr=${{#chr[@]}}
        if [[ ${{n_chr}} -eq 1 ]]
        then
            cp {input} {output}
        else
            bcftools concat -Oz -o {output} {input}
        fi

        tabix -f {output}
        '''


#--------------CALCULATE-----------#
#--------------ACCURACY------------#


rule NGSRELATE_IMPUTED:
    input:
        os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz")
    output:
        RESULTS = os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "IMPUTED", config['PREFIX_NAME'] + "_" + str.capitalize("{tool}") + ".ngsrelate")
    params:
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "IMPUTED")
    threads: workflow.cores * 0.5
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "NGSRELATE_IMPUTED_{tool}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1

        id_list=$(mktemp -p {params.OUTDIR})
        bcftools query -l {input} > $id_list
        ngsRelate -h {input} -T GT -z $id_list -O {output.RESULTS} -p {threads}
        rm $id_list
        '''

rule KGD_IMPUTED:
    input:
        os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz")
    output:
        RELATEDNESS = os.path.join(config['OUTPUT_DIR'], "KGD", "IMPUTED", "{tool}", config['PREFIX_NAME'] + "_{tool}_pairwise_relatedness.txt"),
        INBREEDING = os.path.join(config['OUTPUT_DIR'], "KGD", "IMPUTED", "{tool}", config['PREFIX_NAME'] + "_{tool}_inbreeding_coefficient.txt")
    params:
        KGD_GMATRIX = os.path.abspath("scripts/KGD/GBS-Chip-Gmatrix.R")
    threads: 1
    script:
        "scripts/KGD.R"


rule GT4HAPPY:
    input:
        TEST = os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz"),
        TRUTH = rules.FILTERVCF_REF.output.TRUTH,
        REF = config['REFERENCE_SEQ']
    output: os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_4compare.vcf.gz")
    params:
        TRUTH_FREQ = os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", "truth"),
        SAMPLE = os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + ".sample"),
        VCF_TMP = os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_TMP.vcf.gz"),
        OUTNAME_HP = config['PREFIX_NAME'] + "_{tool}_accuracy"
    threads: 1
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "COMPARE_{tool}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        # tabix -f {input.TEST}
        bcftools query -l {input.TRUTH} > {params.SAMPLE}
        bcftools annotate -x ^FORMAT/GT {input.TEST} | bcftools view -S {params.SAMPLE} -Oz -o {params.VCF_TMP}
        tabix -f {params.VCF_TMP}
        bcftools view -R {input.TRUTH} -Oz -o {output} {params.VCF_TMP}
        tabix -f {output}
        rm {params.VCF_TMP}*
        vcftools --gzvcf {input.TRUTH} --freq2 --out {params.TRUTH_FREQ}
        '''

rule VARIANT_ACCURACY:
    input:
        TEST = os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz"),
        TRUTH = rules.FILTERVCF_REF.output.TRUTH
    output:
        multiext(os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME']), "_{tool}_accuracy_metrics_per_variant.tsv", "_{tool}_accuracy_metrics_per_sample.tsv")
        # SAMPLE = os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy_metrics_per_sample.tsv"),
        # VARIANT = os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy_metrics_per_variant.tsv")
    threads: workflow.cores * 0.5
    script:
        "scripts/accuracy_metrics.R"


rule HAPLOTYPE_ACCURACY:
    input:
        TEST = rules.GT4HAPPY.output,
        TRUTH = rules.FILTERVCF_REF.output.TRUTH,
        REF = config['REFERENCE_SEQ']
    output: os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy.extended.csv")
    params:
        OUTNAME = config['PREFIX_NAME'] + "_{tool}_accuracy"
    threads: workflow.cores * 0.5
    container: "docker://pkrusche/hap.py"
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "HAPLOTYPE_ACCURACY_{tool}.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        cd $(dirname {output})
        /opt/hap.py/bin/hap.py {input.TRUTH} {input.TEST} --threads {threads} -o {params.OUTNAME} -r {input.REF}
        '''


## read pop file
popfile = pd.read_csv(config['STAIRWAY_PLOT']['POPFILE'], sep="\t")

# Count unique values in the second column
NPOP = popfile.iloc[:, 1].nunique()
POPS = []

for i in range(NPOP):  # Fix loop function
    POPS.append(str(popfile.iloc[:, 1].unique()[i]))  # Append to a list

POPS = " ".join(POPS) 

rule FST:
    input: 
        IMPUTE = expand(os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz"), tool = IMPUTATION),
        TRUTH = rules.FILTERVCF_REF.output.TRUTH,
        POPFILE = config['STAIRWAY_PLOT']['POPFILE']
    output: os.path.join(config['OUTPUT_DIR'], "POPGEN", "FST", config['PREFIX_NAME'] + "_FST.txt")
    params: 
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "POPGEN", "FST")
    threads: 1
    log: os.path.join(config['OUTPUT_DIR'], "LOG", "FST.log")
    shell:
        '''
        exec > >(tee "{log}") 2>&1
        fst_param=""
        for s in {POPS}; do
            SUBFILE={params.OUTDIR}/subpop${{s}}.txt
            awk -F'\t' -v val=$s '$2 == val' {input.POPFILE} > $SUBFILE
            fst_param="${{fst_param}}--weir-fst-pop ${{SUBFILE}} "
        done

        echo -e "tImputation\tPairwise\tFst" > fst.tsv
        cd {params.OUTDIR}
        for s1 in {POPS}; do
            for s2 in {POPS}; do
                if [[ ! $s1 < $s2 ]]; then
                    continue
                else
                    vcftools --gzvcf {input.TRUTH} \
                    --weir-fst-pop subpop${{s1}}.txt \
                    --weir-fst-pop subpop${{s2}}.txt \
                    --out TRUTH_subpop_${{s1}}_${{s2}} > TRUTH_subpop_${{s1}}_${{s2}}.log 2>&1
                    fst=$(grep "mean" TRUTH_subpop_${{s1}}_${{s2}}.log | cut -d" " -f7)
                    echo -e "TRUTH\t${{s1}}-${{s2}}\t${{fst}}" >> {output}
                    
                    for i in {input.IMPUTE}; do
                        tool=$(basename $(dirname $i))
                        vcftools --gzvcf $i \
                        --weir-fst-pop subpop${{s1}}.txt \
                        --weir-fst-pop subpop${{s2}}.txt \
                        --out ${{tool}}_subpop_${{s1}}_${{s2}} > ${{tool}}_subpop_${{s1}}_${{s2}}.log 2>&1
                        fst=$(grep "mean" ${{tool}}_subpop_${{s1}}_${{s2}}.log | cut -d" " -f7)
                        echo -e "${{tool}}\t${{s1}}-${{s2}}\t${{fst}}" >> {output}
                    done
                fi
            done
        done
        # rm -f {output}
        # for vcf in {input.IMPUTE}; do
        #     echo ${{vcf}}
        #     tool=$(basename $(dirname ${{vcf}}))
        #     echo ${{tool}}
        #     vcftools --gzvcf ${{vcf}} ${{fst_param}} --out ${{tool}} > ${{tool}}.log 2>&1
        #     fst=$(grep "mean" ${{tool}}.log | cut -d" " -f7)
        #     echo -e "${{tool}}\t${{fst}}" >> {output}
        # done
        
        # vcftools --gzvcf {input.TRUTH} ${{fst_param}} --out TRUTH > TRUTH.log 2>&1
        # fst=$(grep "mean" TRUTH.log | cut -d" " -f7)
        # echo -e "TRUTH\t${{fst}}" >> {output}
        
        '''



rule STAIRWAY_IMPUTED:
    input: 
        VCF = os.path.join(config['OUTPUT_DIR'], "{tool}", config['PREFIX_NAME'] + "_{tool}.vcf.gz"),
        POPFILE = config['STAIRWAY_PLOT']['POPFILE']
    output:
        os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "{tool}", "stairway_results.txt")
        # expand(os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "{tool}", "{subpop}", config['PREFIX_NAME'] + "_{tool}_{subpop}_STAIRWAY.png"), subpop = POPS, tool = "{tool}")
    params:
        BLUEPRINT = "scripts/template.blueprint",
        SEQ_LEN = str(config['STAIRWAY_PLOT']['SEQUENCE_LENGTH']),
        MU = config['STAIRWAY_PLOT']['MUT_RATE'],
        PLOT_TITLE = "{tool}-" + config['PREFIX_NAME'],
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "{tool}"),
        PREFIX_NAME = config['PREFIX_NAME']
    script: "scripts/run_SFS_Stairway.R"


rule STAIRWAY_TRUTH:
    input: 
        VCF = rules.FILTERVCF_REF.output.TRUTH,
        POPFILE = config['STAIRWAY_PLOT']['POPFILE']
    output:
        os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "TRUTH", "stairway_results.txt")
        # expand(os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "{tool}", "{subpop}", config['PREFIX_NAME'] + "_{tool}_{subpop}_STAIRWAY.png"), subpop = POPS, tool = "{tool}")
    params:
        BLUEPRINT = "scripts/template.blueprint",
        SEQ_LEN = config['STAIRWAY_PLOT']['SEQUENCE_LENGTH'],
        MU = config['STAIRWAY_PLOT']['MUT_RATE'],
        PLOT_TITLE = "Groundtruth-"+ config['PREFIX_NAME'],
        OUTDIR = os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "TRUTH"),
        PREFIX_NAME = config['PREFIX_NAME']
    script: "scripts/run_SFS_Stairway.R"


rule POPGEN:
    input: 
        expand(rules.STAIRWAY_IMPUTED.output, tool = IMPUTATION),
        os.path.join(config['OUTPUT_DIR'], "POPGEN", "STAIRWAY", "TRUTH", "stairway_results.txt"),
        rules.FST.output if NPOP > 1 else ()


rule ACCURACY:
    input:
        NGSRELATE_BAM = rules.NGSRELATE_BAM.output,
        NGSRELATE_TRUTH = rules.NGSRELATE_TRUTH.output,
        NGSRELATE_IMPUTED = expand(os.path.join(config['OUTPUT_DIR'], "NGSRELATE", "IMPUTED", config['PREFIX_NAME'] + "_{tool}.ngsrelate"), tool = IMPUTATION),
        # KGD_LCS = rules.KGD_LCS.output.RELATEDNESS,
        # KGD_TRUTH = rules.KGD_TRUTH.output.RELATEDNESS,
        # KGD_IMPUTED = expand(os.path.join(config['OUTPUT_DIR'], "KGD", "IMPUTED", "{tool}", config['PREFIX_NAME'] + "_{tool}_pairwise_relatedness.txt"), tool = IMPUTATION),
        METRICS_SAM = expand(os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy_metrics_per_sample.tsv"), tool = IMPUTATION),
        METRICS_VAR = expand(os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy_metrics_per_variant.tsv"), tool = IMPUTATION),
        HAPPY = expand(os.path.join(config['OUTPUT_DIR'], "ACCURACY", "{tool}", config['PREFIX_NAME'] + "_{tool}_accuracy.extended.csv"), tool = IMPUTATION),
        POPGEN = rules.POPGEN.input
    output:
        NGSRELATE_REL = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_Relatedness.png"),
                # category="Output (Imputed data)",
                labels={"figure": "1 - ngsRelate Pairwise Relatedness"}
            ),
        NGSRELATE_INBREEDING = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_Inbreeding.png"),
                # category="Output (Imputed data)",
                labels={"figure": "2 - ngsRelate Inbreeding Coefficient"}
            ),
        # KGD = report(
        #         os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_G5.png"),
        #         labels={"figure": "3 - KGD Relatedness"}
        #     ),
        HAPPY = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_happy.png"),
                labels={"figure": "6 - Hap.py Results"}
            ),
        VARIANT = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_accuracy_per_variant.png"),
                labels={"figure": "7 - Accuracy Per Variant"}
            ),
        SAMPLE = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_accuracy_per_sample.png"),
                labels={"figure": "8 - Accuracy Per Sample"}
            ),
        MAF = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_MAF_dist.png"),
                labels={"figure": "9 - Truth vs. Imputed MAF"}
            ),
        BENCHMARK = report(
                os.path.join(config['OUTPUT_DIR'], "ACCURACY", "PLOT", config['PREFIX_NAME'] + "_benchmark.png"),
                labels={"figure": "10 - Benchmark"}
            )
    params:
        BENCHMARK = os.path.join(config['OUTPUT_DIR'], "BENCHMARK")
    threads: 1
    script:
        "scripts/plot_imputation_accuracy.R"

