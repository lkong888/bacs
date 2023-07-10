"""
Workflow to compute taps conversion(notrim) #####merge bam file of both 26May2023 and 23June2023
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 8 --snakefile code/bacs_mRNA_merge_spikein.smk --cluster "sbatch -p short -o bacs2.log -e bacs2.err --job-name=bacs2 --constraint=hga --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["merge1"]
print(SAMPLES)
nn_fa = config["nn_fa"]
nn_fasta = config["nn_fasta"]

rule all:
    input:
        expand("mergebam/{sample}.nnunn.clean.filter.sort_pseusite.mpile.txt", sample = SAMPLES),
        expand("mergebam/{sample}.rRNA.filter.dedup_pseusite.mpile.txt", sample = SAMPLES),
        expand("mergebam/{sample}.clean.nnunn.sort_NNUNN1.contex.table", sample = SAMPLES)

#######################################################################merge spike in and nnunn #############################################################
rule merge_nnunn: 
    input:
        bam1=expand("/users/ludwig/ebu571/ebu571/26May2023/align/{sample}.clean.nnunn.sort.bam",  sample=['H1_2023May26_S1','H2_2023May26_S2', 'H3_2023May26_S3', 'H4_2023May26_S4','H5_2023May26_S5', 'H6_2023May26_S6','H7_2023May26_S7','H8_2023May26_S8']),
        bam2=expand("align/{sample}.clean.nnunn.sort.bam",  sample=['H1_2023Jun19_S1','H2_2023Jun19_S2', 'H3_2023Jun19_S3', 'H4_2023Jun19_S4','H5_2023Jun19_S5', 'H6_2023Jun19_S6','H7_2023Jun19_S7','H8_2023Jun19_S8'])
    output:
        expand("mergebam/{sample}.clean.nnunn.sort.bam",sample=['H1','H2', 'H3','H4','H5','H6', 'H7','H8']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[0]} {input.bam1[0]} {input.bam2[0]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[1]} {input.bam1[1]} {input.bam2[1]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[2]} {input.bam1[2]} {input.bam2[2]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[3]} {input.bam1[3]} {input.bam2[3]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[4]} {input.bam1[4]} {input.bam2[4]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[5]} {input.bam1[5]} {input.bam2[5]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[6]} {input.bam1[6]} {input.bam2[6]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[7]} {input.bam1[7]} {input.bam2[7]}
        """
rule nnunn_count:
    input:
        "mergebam/{sample}.clean.nnunn.sort.bam"
    output:
        mpile="mergebam/{sample}.nnunn.clean.filter.sort.mpile.txt",
        sta="mergebam/{sample}.nnunn.clean.filter.sort_pseusite.mpile.txt"
    params:
        ref=nn_fasta
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {input}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {input} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15}}' > {output.sta}
        """

rule nnunn_context:
    input:
        "mergebam/{sample}.clean.nnunn.sort.bam"
    output:
        "mergebam/{sample}.clean.nnunn.sort_NNUNN1.contex.table"
    envmodules: "R/4.0.3-foss-2020b "
    shell:
        """
        Rscript code/nnunn_sta.r {input}
        """


#############################merge spikein######################################################
rule merge_spikein: 
    input:
        bam1=expand("/users/ludwig/ebu571/ebu571/26May2023/align/{sample}.rRNA.filter.dedup.bam",  sample=['H1_2023May26_S1','H2_2023May26_S2', 'H3_2023May26_S3', 'H4_2023May26_S4','H5_2023May26_S5', 'H6_2023May26_S6','H7_2023May26_S7','H8_2023May26_S8']),
        bam2=expand("align/{sample}.rRNA.filter.dedup.bam",  sample=['H1_2023Jun19_S1','H2_2023Jun19_S2', 'H3_2023Jun19_S3', 'H4_2023Jun19_S4','H5_2023Jun19_S5', 'H6_2023Jun19_S6','H7_2023Jun19_S7','H8_2023Jun19_S8'])
    output:
        expand("mergebam/{sample}.rRNA.filter.dedup.bam",sample=['H1','H2', 'H3','H4','H5','H6', 'H7','H8']),
        
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[0]} {input.bam1[0]} {input.bam2[0]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[1]} {input.bam1[1]} {input.bam2[1]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[2]} {input.bam1[2]} {input.bam2[2]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[3]} {input.bam1[3]} {input.bam2[3]} 
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[4]} {input.bam1[4]} {input.bam2[4]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[5]} {input.bam1[5]} {input.bam2[5]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[6]} {input.bam1[6]} {input.bam2[6]}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output[7]} {input.bam1[7]} {input.bam2[7]}
        """

rule merge_count: 
    input:
        "mergebam/{sample}.rRNA.filter.dedup.bam"
    output:
        mpile="mergebam/{sample}.rRNA.filter.dedup.mpile.txt",
        sta="mergebam/{sample}.rRNA.filter.dedup_pseusite.mpile.txt"
    params:
        ref="resource/spikein_rRNA.fa"
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {input} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | grep ^spike_in | awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2,  $3, $4, $5+$14, $7+$16, $11+$20}}' > {output.sta}
        cat {output.mpile} | grep -v ^spike_in | grep -v ^nnunn| awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' >> {output.sta}
        
        """

