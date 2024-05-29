###This script shows an example of calling pseudouridine sites
"""
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 14 --snakefile code/bacs_smallrna_callsite.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs --constraint=hga --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

sample = config["ribo_depleted"]
print(sample)
treat = ['H13_2024May09_S13', 'H14_2024May09_S14']
ctrl = ['H21_2024May09_S21', 'H22_2024May09_S22']

INDEX = ["13","14"]

rule all:
    input:
        expand("smallrna/treat{index}_input{index}_rrna_smallrna_allsites.txt", index=INDEX),
        expand("align/{sample}.rRNA.snoRNA.tRNA.filter.sort_pseusite.mpile.txt", sample=sample)


##############################################################################################################################################################################
rule combine_site:
  input:
    rrna="align/{sample}.rRNA.filter.dedup_pseusite.mpile.txt",
    trna="align/{sample}.tRNA.filter.sort_pseusite.mpile.txt",
    snorna="align/{sample}.snoRNA.filter.sort_pseusite.mpile.txt"
  output:
    "align/{sample}.rRNA.snoRNA.tRNA.filter.sort_pseusite.mpile.txt"
  threads: 8
  envmodules:
    "BEDTools/2.30.0-GCC-10.2.0"
  shell:
     """
    cat {input.rrna} {input.trna} {input.snorna} > {output}
    rm {input.rrna}
    rm {input.snorna}
    rm {input.trna}
     """

rule conversion_fp:
  input:
    treat_nnunn1=expand("align/{sample}.clean.nnunn.sort_NNUNN1.contex.table", sample=treat),
    treat_nnunn2=expand("align/{sample}.clean.nnunn.sort_NNUNN2.contex.table", sample=treat),
    ctrl_nnunn1=expand("align/{sample}.clean.nnunn.sort_NNUNN1.contex.table", sample=ctrl),
    ctrl_nnunn2=expand("align/{sample}.clean.nnunn.sort_NNUNN2.contex.table", sample=ctrl)
  output:
    expand("smallrna/treat{index}_input{index}_spikein.table",index=INDEX)
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/conversion_fp.r {input.ctrl_nnunn1[0]} {input.treat_nnunn1[0]} {input.ctrl_nnunn2[0]} {input.treat_nnunn2[0]} {output[0]}
        Rscript code/conversion_fp.r {input.ctrl_nnunn1[1]} {input.treat_nnunn1[1]} {input.ctrl_nnunn2[1]} {input.treat_nnunn2[1]} {output[1]}
     """

rule group:
  input:
    treat=expand("align/{sample}.rRNA.snoRNA.tRNA.filter.sort_pseusite.mpile.txt", sample=treat),
    ctrl=expand("align/{sample}.rRNA.snoRNA.tRNA.filter.sort_pseusite.mpile.txt", sample=ctrl)
  params:
    bed=expand("align/treat{index}.bed",index=INDEX),
    fabed=expand("align/treat{index}.fasta.bed",index=INDEX),
    sta=expand("align/treat{index}.sta.txt",index=INDEX),
    ref="/users/ludwig/ebu571/ebu571/project/psi_ko/resource/rRNA_tRNA_snoRNA.fa"
  output:
    groupbed=expand("smallrna/treat{index}_input{index}_rrna_smallrna_allT_motif.bed",index=INDEX)
  threads: 8
  envmodules:
    "BEDTools/2.30.0-GCC-10.2.0"
  shell:
     """
        cat {input.treat[0]} | awk '$2 > 3' | awk '{{OFS="\\t"}}{{print $1, $2-3, $2+2, $3, $4, $5, $6, $7,$8}}' > {params.bed[0]}
        bedtools getfasta -fi {params.ref} -bed {params.bed[0]} | paste - - | sed 's/>//g' | sed 's/:/\\t/g' | sed 's/-/\\t/g' | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4}}' > {params.fabed[0]}
        cat {params.bed[0]} | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4, $5, $6, $7, $8,$9}}' | bedtools intersect -a - -b {params.fabed[0]} -wa -wb -loj | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8,$9, $13}}' > {params.sta[0]}
        bedtools intersect -a <(cat {params.sta[0]} | sort -k1,1 -k2,2n) -b <(cat {input.ctrl[0]} | sed 's/-/_/g' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, $3, $4, $5, $6, $7,$8}}' | sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($6 + $7) > 0' | awk '($16 + $17) > 0' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $7/($6 + $7), $8, $8/$5, $9,$10, $15, $16,$17, $17/($16 + $17), $18, $18/$15,$19}}' > {output.groupbed[0]}

        cat {input.treat[1]} | awk '$2 > 3' | awk '{{OFS="\\t"}}{{print $1, $2-3, $2+2, $3, $4, $5, $6, $7,$8}}' > {params.bed[1]}
        bedtools getfasta -fi {params.ref} -bed {params.bed[1]} | paste - - | sed 's/>//g' | sed 's/:/\\t/g' | sed 's/-/\\t/g' | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4}}' > {params.fabed[1]}
        cat {params.bed[1]} | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4, $5, $6, $7, $8,$9}}' | bedtools intersect -a - -b {params.fabed[1]} -wa -wb -loj | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8,$9, $13}}' > {params.sta[1]}
        bedtools intersect -a <(cat {params.sta[1]} | sort -k1,1 -k2,2n) -b <(cat {input.ctrl[1]} | sed 's/-/_/g' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, $3, $4, $5, $6, $7,$8}}' | sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($6 + $7) > 0' | awk '($16 + $17) > 0' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $7/($6 + $7), $8, $8/$5, $9,$10, $15, $16,$17, $17/($16 + $17), $18, $18/$15,$19}}' > {output.groupbed[1]}

        rm {params.bed}
        rm {params.fabed}
        rm {params.sta}

     """

rule test:
  input:
    spikein="smallrna/treat{index}_input{index}_spikein.table",
    groupbed="smallrna/treat{index}_input{index}_rrna_smallrna_allT_motif.bed",

  output:
    "smallrna/treat{index}_input{index}_rrna_smallrna_allsites.txt",

  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/calling_rRNA_snorna_trna_update.r {input.spikein} {input.groupbed} {output}
     """

