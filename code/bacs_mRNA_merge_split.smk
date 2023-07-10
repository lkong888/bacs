"""
Workflow to compute taps conversion(notrim) #####merge bam file of both 26May2023 and 23June2023
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 8 --snakefile code/bacs_mRNA_merge1.smk --cluster "sbatch -p short -o bacs_split.log -e bacs_split.err --job-name=bacssplit --constraint=hga --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["merge"]
print(SAMPLES)
CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"] 

rule all:
    input:
        #expand("mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam", sample = SAMPLES),
        #expand("mergebam/{sample}.mRNA_first.mpile.txt", sample = SAMPLES),
        #expand("mergebam/{sample}.mapping_report.txt", sample = SAMPLES),
        expand("splitbam/{sample}.mRNA_pseusite.mpile.txt.gz", sample = SAMPLES),
        #("rpkm/{sample}-chrM.tpm", sample = SAMPLES)

rule merge: 
    input:
        bam1=expand("/users/ludwig/ebu571/ebu571/26May2023/align/{sample}.mRNAAligned.sortedByCoord.out.bam",  sample=['H1_2023May26_S1','H2_2023May26_S2', 'H3_2023May26_S3', 'H4_2023May26_S4','H5_2023May26_S5', 'H6_2023May26_S6','H7_2023May26_S7','H8_2023May26_S8']),
        bam2=expand("align/{sample}.mRNAAligned.sortedByCoord.out.bam",  sample=['H1_2023Jun19_S1','H2_2023Jun19_S2', 'H3_2023Jun19_S3', 'H4_2023Jun19_S4','H5_2023Jun19_S5', 'H6_2023Jun19_S6','H7_2023Jun19_S7','H8_2023Jun19_S8'])
    output:
        expand("mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam",sample=['treat1','treat2', 'input1','input2']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[0]} {input.bam1[0]} {input.bam2[0]} {input.bam1[1]} {input.bam2[1]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[1]} {input.bam1[2]} {input.bam2[2]} {input.bam1[3]} {input.bam2[3]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[2]} {input.bam1[4]} {input.bam2[4]} {input.bam1[5]} {input.bam2[5]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[3]} {input.bam1[6]} {input.bam2[6]} {input.bam1[7]} {input.bam2[7]}

        """

rule collect_metrics:
    input:
        "mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam"
    output:
        "mergebam/{sample}.star_mRNA.collect_metrics.alignment_summary_metrics"
    params:
        prefix="mergebam/{sample}.star_mRNA.collect_metrics",
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
        "mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam"
    output:
        sortbam="mergebam/{sample}.mrna.filter.sort.bam",
        dedup="mergebam/{sample}.mrna.filter_dedup.sort.bam"
    params:
        mq=" -q 30",
        tmp="{sample}.tmp",
        prefix="logs/{sample}.mRNA.dedup"
    envmodules: "UMI-tools/1.0.1-foss-2019b-Python-3.7.4"
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {params.mq} {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=10' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}
        umi_tools dedup --random-seed=123 --method=unique -I {output.sortbam} -S {output.dedup} --output-stats {params.prefix}.stats -L {params.prefix}.logs
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.dedup}
        """


rule split: 
    input:
        "mergebam/{sample}.mrna.filter_dedup.sort.bam"
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

rule merge1: 
    input:
        first="splitbam/{sample}.mRNA.{chrs}_first.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.bam"
    output:
        first="splitbam/{sample}.mRNA.{chrs}_first.mpile.txt.gz",
        second="splitbam/{sample}.mRNA.{chrs}_second.mpile.txt.gz",
        merge="splitbam/{sample}.mRNA.{chrs}_pseusite.mpile.txt.gz",

    envmodules:  "CMake/3.23.1-GCCcore-11.3.0"
    threads: 2
    shell:
        """
            /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -AB -d 0 -Q 10 --reverse-del \
            -f /users/ludwig/ebu571/ebu571/resource/GRCh38.fa {input.first} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | gzip > {output.first}

            /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -AB -d 0 -Q 10 --reverse-del \
            -f /users/ludwig/ebu571/ebu571/resource/GRCh38.fa {input.second} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | gzip > {output.second}
            
            zcat {output.first}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>5'| awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "-", $3, $4, $5+$14, $7+$16}}' | gzip > {output.merge}
            zcat {output.second}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>5' | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "+", $3, $4, $8+$17, $6+$15}}' | gzip >> {output.merge}
            
        """

rule sta: 
    input:
        first1=expand("splitbam/{{sample}}.mRNA.{chrs}_first.mpile.txt.gz", chrs=CHRS), 
        second1=expand("splitbam/{{sample}}.mRNA.{chrs}_second.mpile.txt.gz", chrs=CHRS)
    output:
        merge1="splitbam/{sample}.mRNA_pseusite.mpile.txt.gz"

    envmodules: "samtools/1.8-gcc5.4.0", "CMake/3.23.1-GCCcore-11.3.0"
    threads: 2
    shell:
        """
            zcat {input.first1}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>5'| awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "-", $3, $4, $5+$14, $7+$16}}' | gzip > {output.merge1}
            zcat {input.second1}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>5' | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "+", $3, $4, $8+$17, $6+$15}}' | gzip >> {output.merge1}
            
        """
