"""
Workflow for merging replicates, featurecounts and mutation calling.
inputs: 
    bam
outputs:
    counts
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 -j 14 --snakefile code/bacs_alignment2.smk --cluster "sbatch -p short -o bacs_merge.log -e bacs_merge.err --job-name=bacs_merge --constraint=hga --cpus-per-task 3 "
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
        expand("align/{merge}.snoRNA.tRNA.filter.sort_A.mpile.txt", merge = INPUT)

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

rule countU:
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
rule countA:
    input:
        "align/{merge}.snoRNA.tRNA.filter.sort.mpile.txt"
    output:
        sta="align/{merge}.snoRNA.tRNA.filter.sort_A.mpile.txt"
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        cat {input} | awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $4-($5+$14)}}' > {output.sta}
        """

