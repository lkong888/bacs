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
        expand("featureCounts/{merge}-chrM.snorna.trna.tpm", merge = INPUT),
        expand("featureCounts/{merge}.tpm", merge = INPUT),
        expand("align/{merge}.snoRNA.tRNA.filter.sort_A.mpile.txt", merge = INPUT),
        "site/snoRNA_called_sites_fdrpassed_metagene.pdf"

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
##cat tRNA_snorna.fa.fai | awk '{printf $1"\t""bacs""\t""tRNA_snorna""\t""1""\t"$2"\t"".""\t"".""\t"".""\t""gene_id ""\""$1"\"""; ""transcript_id ""\""$1"\""";""\n"}' > snorna_trna.gtf
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


rule combine_counts_rpkms:
  input:
    counts = expand("featureCounts/{merge}.snorna.trna.counts",merge = INPUT),
    tpms = expand("featureCounts/{merge}-chrM.snorna.trna.tpm", merge = INPUT)
  output: 
    counts = "featureCounts/{merge}.counts",
    tpm = "featureCounts/{merge}.tpm",
  threads: 1
  shell:
    "len=$(ls {input.counts}|wc -l);"
    "paste {input.tpms} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) > {output.tpm};"
    "paste {input.counts} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) |"
    "grep -v 'chrM' > {output.counts}"


rule deeptools:
    input:
        "site/snoRNA_tRNA_called_sites_fdrpassed.txt"
    output:
        bg="site/snoRNA_called_sites_fdrpassed.bedGraph",
        bw="site/snoRNA_called_sites_fdrpassed.bw",
        mat="site/snoRNA_called_sites_fdrpassed.mat.gz",
        pdf="site/snoRNA_called_sites_fdrpassed_metagene.pdf"
    envmodules: 
        "deepTools/3.3.1-foss-2018b-Python-3.6.6", "plotly.py/4.4.1-foss-2018b-Python-3.6.6"
    params:
        region="resource/snoRNA_reference.fasta.bed",
        gsize="resource/snoRNA_reference.size.txt"
    shell: 
        """
        cat {input} | grep -v tRNA | grep -v ^chr | awk '{{OFS="\\t"}}{{print $1,$2,$3,($8+$16+$22)/3}}' | sort -k1,1 -k2,2n > {output.bg}

        /well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig {output.bg} {params.gsize} {output.bw}

        computeMatrix scale-regions -S {output.bw} -R {params.region} \
                              --regionBodyLength 100 \
                              -o {output.mat} --binSize 10 -p 8

        plotProfile -m {output.mat} \
              -out {output.pdf} 
        
        """

#####cat snoRNA_reference.fasta.fai | awk '{OFS="\t"}{print $1, $2}'> snoRNA_reference.size.txt
#####cat snoRNA_reference.fasta.fai | awk '{OFS="\t"}{print $1, "0", $2}' > snoRNA_reference.fasta.bed