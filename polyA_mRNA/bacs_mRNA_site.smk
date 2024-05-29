"""
Workflow to compute taps conversion(notrim) #####merge bam file of both 26May2023 and 23June2023
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 8 --snakefile code/bacs_mRNA_site.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs1 --cpus-per-task 8 "
"""
rule all:
    input:
        expand("site/treat_input_allT_motif{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("site/treat_input_pvalue{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("site/treat_input_ivt_pvalue{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("site/treat_input_ivt_pvalue1{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09']),

##############call site
rule motif:
  input:
    treat="splitbam/treat.mRNA_pseusite.mpile.txt.gz",
  output:
    sta="site/treat.sta.txt",
  params:
    bed="site/treat.bed",
    fabed="site/treat.fasta.bed",
  threads: 8
  envmodules:
    "BEDTools/2.30.0-GCC-10.2.0"
  shell:
     """
        zcat {input.treat}| awk '$2>3 && ($7+$8)>=20' |awk '{{OFS="\\t"}}{{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}}' > {params.bed}
        bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed {params.bed} | paste - - | sed 's/>//g' | sed 's/:/\\t/g'| sed 's/-/\\t/g' | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4}}'> {params.fabed}
        zcat {input.treat} | bedtools intersect -a - -b {params.fabed} -wa -wb -loj | awk '($7+$8)>=20' |awk '{{OFS="\\t"}}{{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}}' > {output.sta}
        rm {params.bed}
        rm {params.fabed}
 """
###to speed-up,only get motif sequence with (T+C)>=20
rule group:
  input:
    treat="splitbam/treat.mRNA_pseusite.mpile.txt.gz",
    ctrl="splitbam/input.mRNA_pseusite.mpile.txt.gz",
    sta="site/treat.sta.txt"
  output:
    groupbed="site/treat_input_allT_motif.bed"
  threads: 8
  envmodules:
    "BEDTools/2.30.0-GCC-10.2.0"
  shell:
     """
        bedtools intersect -a <(cat {input.sta} |sort -k1,1 -k2,2n) -b <(zcat {input.ctrl}| sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>=20 && ($17+$18)>=20'| awk '{{OFS="\\t"}}{{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$9/$6,$10,$15,$16,$17,$18,$18/($17+$18),$19,$19/$16}}' > {output.groupbed}
     """
###only (T+C)>=20 in both treat and control will be output
rule splitbed:
  input:
    "site/treat_input_allT_motif.bed"
  output:
    expand("site/treat_input_allT_motif{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09'])
  params:
    prefix="site/treat_input_allT_motif"
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
     """
        linenum0=`echo "($(wc -l < "{input}") / 10) + 1" | bc`
        split -d {input} {params.prefix} -l $linenum0
        
     """

rule test:
  input:
    expand("site/treat_input_allT_motif{{splitid}}")
  output:
    expand("site/treat_input_pvalue{{splitid}}")
  params:
    spikein="mergebam/mrna_spikein.table"
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/calling_mRNA.r {input} {output} {params.spikein}
     """

rule test_ivt1:
  input:
    expand("site/treat_input_allT_motif{{splitid}}"),
  output:
    expand("site/treat_input_ivt_pvalue1{{splitid}}")
  params:
    spikein="mergebam/mrna_spikein.table"
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/calling_mRNA_ivt1.r {input} {output} splitbam/ivt.mRNA_pseusite.mpile.abover10.txt {params.spikein}
     """
  
