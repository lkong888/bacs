"""
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 -j 8 --snakefile code/bacs_mRNA_part2.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs1 --cpus-per-task 8 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["mrna"]
print(SAMPLES)
#CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"] 
CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"] 

rule all:
    input:
        expand("align/{sample}.mRNAAligned.sortedByCoord.out.bam", sample = SAMPLES),
        expand("align/{sample}.mapping_report.txt", sample = SAMPLES),
        "rpkm/mrna.tpm",
        expand("splitbam/{sample}.mRNA_pseusite.mpile.txt.gz", sample = SAMPLES),

#######################################################################merge mRNA#############################################################
rule merge: 
    input:
        bam=expand("align/{sample}.mRNAAligned.sortedByCoord.out.bam",  sample=['H1_2023Oct26_S1','H2_2023Oct26_S2', 'H3_2023Oct26_S3', 'H4_2023Oct26_S4', 'H5_2023Oct26_S5', 'H6_2023Oct26_S6','H7_2023Oct05_S7','H8_2023Oct05_S8'])
    output:
        expand("align/{sample}.mRNAAligned.sortedByCoord.out.bam",sample=['input','treat', 'ivt']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[0]} {input.bam[0]} {input.bam[1]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[1]} {input.bam[2]} {input.bam[3]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[2]} {input.bam[4]} {input.bam[5]} {input.bam[6]} {input.bam[7]}

        """

rule collect_metrics:
    input:
        "align/{sample}.mRNAAligned.sortedByCoord.out.bam"
    output:
        "align/{sample}.star_mRNA.collect_metrics.alignment_summary_metrics"
    params:
        prefix="align/{sample}.star_mRNA.collect_metrics",
        ref="/users/ludwig/ebu571/ebu571/resource/GRCh38.fa"
    envmodules: 
        "picard/2.23.0-Java-11", "R/4.0.3-foss-2020b"
    shell:
        """
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectMultipleMetrics \
            I={input}  \
            O={params.prefix} \
            R={params.ref}
        """

rule star_dedup: 
    input:
        "align/{sample}.mRNAAligned.sortedByCoord.out.bam"
    output:
        sortbam="align/{sample}.mrna.filter.sort.bam",
        dedup="align/{sample}.mrna.filter_dedup.sort.bam"
    params:
        mq=" -q 30",
        tmp="{sample}.tmp",
        prefix="logs/{sample}.mRNA.dedup"
    envmodules: "UMI-tools/1.0.1-foss-2019b-Python-3.7.4"
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {params.mq} {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=3' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}
        umi_tools dedup --random-seed=123 --method=unique -I {output.sortbam} -S {output.dedup} --output-stats {params.prefix}.stats -L {params.prefix}.logs
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.dedup}
        """

rule mapping:
    input:
        mergebam="align/{sample}.mRNAAligned.sortedByCoord.out.bam",
        filt_dedup="align/{sample}.mrna.filter_dedup.sort.bam"
    output:
        "align/{sample}.mapping_report.txt", 
    envmodules: "samtools/1.8-gcc5.4.0"
    shell:
        """
        echo -e "smp\\tmergebam\\tfilt_dedup" > {output}
        mergebam=`samtools view {input.mergebam} | cut -f 1 | sort | uniq | wc -l`
        filt_dedup=`samtools view {input.filt_dedup} | cut -f 1 | sort | uniq | wc -l`
        echo -e "{output}\\t${{mergebam}}\\t${{filt_dedup}}" >> {output}
        """
rule feature_count_rpkm:
  input:
    bam ="align/{sample}.mrna.filter_dedup.sort.bam",
  output:
    counts = "featureCounts/{sample}.counts",
    tpm = "rpkm/{sample}-chrM.tpm"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
    "featureCounts -a /gpfs3/well/ludwig/users/ebu571/resource/gencode.v43.annotation.gtf -o {output.counts} {input.bam} "
    " -F GTF -T {threads} -t exon -g gene_name; "
    "grep -v 'chrM' {output.counts} | Rscript /users/ludwig/ebu571/ebu571/20May2023_final/code/featurecounts2tpm.r - > {output.tpm}"

rule combine_counts_rpkms:
  input:
    counts = expand("featureCounts/{sample}.counts",sample = SAMPLES),
    tpms = expand("rpkm/{sample}-chrM.tpm", sample = SAMPLES)
  output: 
    counts = "rpkm/mrna.counts",
    tpm = "rpkm/mrna.tpm",
  threads: 1
  shell:
    "len=$(ls {input.counts}|wc -l);"
    "paste {input.tpms} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) > {output.tpm};"
    "paste {input.counts} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) |"
    "grep -v 'chrM' > {output.counts}"



rule split: 
    input:
        "align/{sample}.mrna.filter_dedup.sort.bam"
    output:
        split=expand("splitbam/{{sample}}.mRNA.{chrs}.sort.bam", chrs=CHRS)
    params:
        mq=" -q 30",
        ref="/users/ludwig/ebu571/ebu571/resource/GRCh38.fa",
        prefix="splitbam/{sample}.mRNA"
    envmodules: "samtools/1.8-gcc5.4.0"
    threads: 2
    shell:
        """
        samtools index {input}
        grep -wP "chr[0-9]*" {params.ref}.fai |cut -f1|while read chr  
        do
        samtools view -bS {input} ${{chr}} >{params.prefix}.${{chr}}.bam
        samtools sort {params.prefix}.${{chr}}.bam > {params.prefix}.${{chr}}.sort.bam
	    done

        grep chr[X,Y,M] {params.ref}.fai |cut -f1|while read chr  
        do
        samtools view -bS {input} ${{chr}} >{params.prefix}.${{chr}}.bam
        samtools sort {params.prefix}.${{chr}}.bam > {params.prefix}.${{chr}}.sort.bam
	    done
        """


rule clasify: 
    input:
        "splitbam/{sample}.mRNA.{chrs}.sort.bam"
    output:
        first="splitbam/{sample}.mRNA.{chrs}_first.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.bam"
    params:
        mq=" -q 30"
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -h {input} |awk '$2==0 ||$1~/^@/' |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - >{output.first}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -h {input} |awk '$2==16 ||$1~/^@/' |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - >{output.second}
        """

#####using polyG for trimming in both strand
rule count:
    input:
        first="splitbam/{sample}.mRNA.{chrs}_first.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.bam"
    output:
        first="splitbam/{sample}.mRNA.{chrs}_first.mpile.txt.gz",
        second="splitbam/{sample}.mRNA.{chrs}_second.mpile.txt.gz",
    envmodules: "CMake/3.23.1-GCCcore-11.3.0"
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -AB -d 0 -Q 10 --reverse-del \
            -f /users/ludwig/ebu571/ebu571/resource/GRCh38.fa {input.first} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | gzip > {output.first}

        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -AB -d 0 -Q 10 --reverse-del \
            -f /users/ludwig/ebu571/ebu571/resource/GRCh38.fa {input.second} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup |gzip > {output.second}
        """


rule sta:
    input:
        first=expand("splitbam/{{sample}}.mRNA.{chrs}_first.mpile.txt.gz", chrs=CHRS), 
        second=expand("splitbam/{{sample}}.mRNA.{chrs}_second.mpile.txt.gz", chrs=CHRS)
    output:
        "splitbam/{sample}.mRNA_pseusite.mpile.txt.gz"
    shell:
        """
        zcat {input.first}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>0'| awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "-", $3, $4, $5+$14, $7+$16,$11+$20}}' | gzip > {output}
        zcat {input.second}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>0' | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "+", $3, $4, $8+$17, $6+$15, $11+$20}}' | gzip >> {output}

        """

