"""
Workflow to compute taps conversion(notrim)
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 -j 14 --snakefile code/analysis_for_bacs_merge.smk --cluster "sbatch -p short -o bacs_merge.log -e bacs_merge.err --job-name=bacs_merge --constraint=hga --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

INPUT = config["group"]
nn_fa = config["nn_fa"]
nn_fasta = config["nn_fasta"]

print(INPUT)

rule all:
    input:
        expand("align/{merge}.snoRNA.tRNA.filter.sort_pseusite.mpile.txt", merge = INPUT),
        expand("align/{merge}.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt", merge = INPUT),
        expand("featureCounts/{merge}-chrM.snorna.trna.tpm", merge = INPUT),
        expand("featureCounts/{merge}-chrM.snorna.trna.uniq.best.tpm", merge = INPUT),
        expand("align/{merge}.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt", merge = INPUT),
        expand("align/{merge}.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt", merge = INPUT),
        expand("featureCounts/{merge}.snorna.trna.mapq10.counts", merge = INPUT),
        expand("featureCounts/{merge}.snorna.trna.mapq1.counts", merge = INPUT)

##############################################################################################################################################################################
############################snorna merge and count ########################

rule merge: 
    input:
        expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam",  sample=['H1_2023May19_S1','H2_2023May19_S2', 'H3_2023May19_S3', 'H4_2023May19_S4','H5_2023May19_S5', 'H6_2023May19_S6','H7_2023May19_S7','H8_2023May19_S8', 'H9_2023May19_S9', 'H10_2023May19_S10','H3_2023Apr16_S3', 'H16_2023Apr16_S15','H23_2023Apr16_S19', 'H26_2023Apr16_S22'])
    output:
        control=expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam",sample=['input1','input2']),
        treat=expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam",sample=['treat1','treat2', 'treat3']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.control[0]} {input[11]} {input[13]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.control[1]} {input[8]} {input[9]}

        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[0]} {input[10]} {input[12]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[1]} {input[0]} {input[2]} {input[4]} {input[6]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[2]} {input[1]} {input[3]} {input[5]} {input[7]}

        """

rule count:
    input:
        "align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam"
    output:
        mpile="align/{merge}.snoRNA.tRNA.filter.sort.mpile.txt",
        sta="align/{merge}.snoRNA.tRNA.filter.sort_pseusite.mpile.txt"
    params:
        ref="resource/tRNA_snorna.fa"
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {input}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {input} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' > {output.sta}
        """
##cat tRNA_snorna.fa.fai | awk '{printf $1"\t""bacs""\t""tRNA_snorna""\t""1""\t"$2"\t"".""\t"".""\t"".""\t""gene_id ""\""$1"\"""; ""transcript_id ""\""$1"\""";""\n"}' > snorna_trna.gtf
rule feature_count:
  input:
    bam ="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam"
  output:
    counts = "featureCounts/{merge}.snorna.trna.counts",
    tpm = "featureCounts/{merge}-chrM.snorna.trna.tpm"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
    "featureCounts -a resource/snorna_trna.gtf -o {output.counts} {input.bam} "
    " -F GTF -T {threads} -t tRNA_snorna -g gene_id; "
    "grep -v 'chrM' {output.counts} | Rscript code/featurecounts2tpm.r - > {output.tpm}"





##########################################results for keeping both uniq reads and best alignmment
rule merge1: 
    input:
        expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam",  sample=['H1_2023May19_S1','H2_2023May19_S2', 'H3_2023May19_S3', 'H4_2023May19_S4','H5_2023May19_S5', 'H6_2023May19_S6','H7_2023May19_S7','H8_2023May19_S8', 'H9_2023May19_S9', 'H10_2023May19_S10','H3_2023Apr16_S3', 'H16_2023Apr16_S15','H23_2023Apr16_S19', 'H26_2023Apr16_S22'])
    output:
        control=expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam",sample=['input1','input2']),
        treat=expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam",sample=['treat1','treat2', 'treat3']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.control[0]} {input[11]} {input[13]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.control[1]} {input[8]} {input[9]}

        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[0]} {input[10]} {input[12]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[1]} {input[0]} {input[2]} {input[4]} {input[6]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.treat[2]} {input[1]} {input[3]} {input[5]} {input[7]}

        """

rule count1:
    input:
        "align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam"
    output:
        mpile="align/{merge}.snoRNA.tRNA.filter.uniq.best.mpile.txt",
        sta="align/{merge}.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt"
    params:
        ref="resource/tRNA_snorna.fa"
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {input}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {input} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' > {output.sta}
        """
##cat tRNA_snorna.fa.fai | awk '{printf $1"\t""bacs""\t""tRNA_snorna""\t""1""\t"$2"\t"".""\t"".""\t"".""\t""gene_id ""\""$1"\"""; ""transcript_id ""\""$1"\""";""\n"}' > snorna_trna.gtf
rule feature_count1:
  input:
    bam ="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam"
  output:
    counts = "featureCounts/{merge}.snorna.trna.uniq.best.counts",
    tpm = "featureCounts/{merge}-chrM.snorna.trna.uniq.best.tpm"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
    "featureCounts -a resource/snorna_trna.gtf -o {output.counts} {input.bam} "
    " -F GTF -T {threads} -t tRNA_snorna -g gene_id; "
    "grep -v 'chrM' {output.counts} | Rscript code/featurecounts2tpm.r - > {output.tpm}"







##########################################################results for MAPQ 10
rule count2:
    input:
        "align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam"
    output:
        filterbam="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.mapq10.bam",
        mpile="align/{merge}.snoRNA.tRNA.filter.mapq10.mpile.txt",
        sta="align/{merge}.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt"
    params:
        ref="resource/tRNA_snorna.fa",
        tmp="{merge}.tmp",
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS -q 10 {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} > {output.filterbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.filterbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {output.filterbam} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' > {output.sta}
        """
##cat tRNA_snorna.fa.fai | awk '{printf $1"\t""bacs""\t""tRNA_snorna""\t""1""\t"$2"\t"".""\t"".""\t"".""\t""gene_id ""\""$1"\"""; ""transcript_id ""\""$1"\""";""\n"}' > snorna_trna.gtf
rule feature_count2:
  input:
    bam ="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.mapq10.bam"
  output:
    counts = "featureCounts/{merge}.snorna.trna.mapq10.counts",
    tpm = "featureCounts/{merge}-chrM.snorna.trna.mapq10.tpm"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
    "featureCounts -a resource/snorna_trna.gtf -o {output.counts} {input.bam} "
    " -F GTF -T {threads} -t tRNA_snorna -g gene_id; "
    "grep -v 'chrM' {output.counts} | Rscript code/featurecounts2tpm.r - > {output.tpm}"


##########################################################results for MAPQ 1
rule count3:
    input:
        "align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam"
    output:
        filterbam="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.mapq1.bam",
        mpile="align/{merge}.snoRNA.tRNA.filter.mapq1.mpile.txt",
        sta="align/{merge}.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt"
    params:
        ref="resource/tRNA_snorna.fa",
        tmp="{merge}.tmp",
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS -q 1 {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} > {output.filterbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.filterbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {output.filterbam} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' > {output.sta}
        """
##cat tRNA_snorna.fa.fai | awk '{printf $1"\t""bacs""\t""tRNA_snorna""\t""1""\t"$2"\t"".""\t"".""\t"".""\t""gene_id ""\""$1"\"""; ""transcript_id ""\""$1"\""";""\n"}' > snorna_trna.gtf
rule feature_count3:
  input:
    bam ="align/{merge}.snoRNA.tRNA.clean.bowtie2.filter.mapq1.bam",
  output:
    counts = "featureCounts/{merge}.snorna.trna.mapq1.counts",
    tpm = "featureCounts/{merge}-chrM.snorna.trna.mapq1.tpm"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
    "featureCounts -a resource/snorna_trna.gtf -o {output.counts} {input.bam} "
    " -F GTF -T {threads} -t tRNA_snorna -g gene_id; "
    "grep -v 'chrM' {output.counts} | Rscript code/featurecounts2tpm.r - > {output.tpm}"
