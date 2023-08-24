"""
Workflow to for merging reads, counting mutation #####merge bam file of both 26May2023 and 23June2023
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 20 --snakefile code/bacs_mRNA_merge.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs1 --cpus-per-task 4 "
"""

from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["merge"]
print(SAMPLES)
CHRS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"] 

rule all:
    input:
        #expand("splitbam/{sample}.mRNA.{chrs}_first.rmsoft1.bam", sample=SAMPLES, chrs=CHRS),
        #expand("splitbam/{sample}.mRNA_pseusite.mpile.txt.gz", sample=SAMPLES),
        #expand("rpkm/{sample}.tpm", sample=SAMPLES)
        #expand("site/treat{index}_input{index}_pvalue{splitid}", index=['1','2'], splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("mergebam/{sample}.mapping_report.txt", sample=SAMPLES)

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
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {params.mq} {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=3' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}
        umi_tools dedup --random-seed=123 --method=unique -I {output.sortbam} -S {output.dedup} --output-stats {params.prefix}.stats -L {params.prefix}.logs
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.dedup}
        """
rule mapping:
    input:
        mergebam="mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam",
        filt_dedup="mergebam/{sample}.mrna.filter_dedup.sort.bam",
    output:
        "mergebam/{sample}.mapping_report.txt", 
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
    bam ="mergebam/{sample}.mrna.filter_dedup.sort.bam",
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
    counts = "rpkm/combine.counts",
    tpm = "rpkm/combine.tpm",
  threads: 1
  shell:
    "len=$(ls {input.counts}|wc -l);"
    "paste {input.tpms} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) > {output.tpm};"
    "paste {input.counts} |cut -f 1-6,$(seq -s, 7 7 $((7*len))) |"
    "grep -v 'chrM' > {output.counts}"



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

rule rmsoft: 
    input:
        first="splitbam/{sample}.mRNA.{chrs}_first.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.bam"
    output:
        rmsoft1="splitbam/{sample}.mRNA.{chrs}_first.rmsoft1.bam",
        rmsoft2="splitbam/{sample}.mRNA.{chrs}_second.rmsoft2.bam",

    log:
        "logs/{sample}.mRNA.{chrs}.rmsoft.log"
    threads: 6
    shell:
        """
        (bash -c '
        . $HOME/.bashrc # if not loaded automatically
        conda activate ngsutils
        bamutils removeclipping {input.first} {output.rmsoft1}
        bamutils removeclipping {input.second} {output.rmsoft2}
        conda deactivate') 2> {log}
       """

rule filterC: 
    input:
        rmsoft1="splitbam/{sample}.mRNA.{chrs}_first.rmsoft1.bam",
        rmsoft2="splitbam/{sample}.mRNA.{chrs}_second.rmsoft2.bam",
    output:
        first="splitbam/{sample}.mRNA.{chrs}_first.clip.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.clip.bam"
    log:
        "logs/{sample}.mRNA.{chrs}.filterC.log"
    envmodules: "GATK/4.1.7.0-GCCcore-8.3.0-Java-11"
    threads: 6
    shell:
        """
        (gatk ClipReads -I {input.rmsoft1} -O {output.first} --disable-tool-default-read-filters true -XF /users/ludwig/ebu571/ebu571/23June2023/resource/polyG.fa -CR WRITE_NS ) 2> {log}
        (gatk ClipReads -I {input.rmsoft2} -O {output.second} --disable-tool-default-read-filters true -XF /users/ludwig/ebu571/ebu571/23June2023/resource/polyG.fa -CR WRITE_NS ) 2>> {log}
        """
rule count:
    input:
        first="splitbam/{sample}.mRNA.{chrs}_first.clip.bam",
        second="splitbam/{sample}.mRNA.{chrs}_second.clip.bam"
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



##############call site
rule motif:
  input:
    treat="splitbam/treat{index}.mRNA_pseusite.mpile.txt.gz"
  output:
    bed="splitbam/treat{index}.bed",
    fabed="splitbam/treat{index}.fasta.bed",
    sta="splitbam/treat{index}.sta.txt",
  threads: 8
  envmodules:
    "BEDTools/2.30.0-GCC-10.2.0"
  shell:
     """
        zcat {input.treat}| awk '$2>3 && ($7+$8)>=20' |awk '{{OFS="\\t"}}{{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}}' > {output.bed}
        bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed {output.bed} | paste - - | sed 's/>//g' | sed 's/:/\\t/g'| sed 's/-/\\t/g' | awk '{{OFS="\\t"}}{{print $1, $2+2, $3-2, $4}}'> {output.fabed}
        zcat {input.treat} | bedtools intersect -a - -b {output.fabed} -wa -wb -loj | awk '($7+$8)>=20' |awk '{{OFS="\\t"}}{{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}}' > {output.sta}
 """
###to speed-up,only get motif sequence with (T+C)>=20
rule group:
  input:
    treat="splitbam/treat{index}.mRNA_pseusite.mpile.txt.gz",
    ctrl="splitbam/input{index}.mRNA_pseusite.mpile.txt.gz",
    sta="splitbam/treat{index}.sta.txt"
  output:
    groupbed="splitbam/treat{index}_input{index}_allT_motif.bed"
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
    expand("splitbam/treat{index}_input{index}_allT_motif.bed",index=["1","2"])
  output:
    expand("splitbam/treat{{index}}_input{{index}}_allT_motif{splitid}", splitid=['00','01','02','03','04','05','06','07','08','09'])
  params:
    prefix=expand("splitbam/treat{index}_input{index}_allT_motif",index=["1","2"])
  threads: 1
  envmodules:
    "Subread/1.6.4-foss-2018b",
    "R/3.6.0-foss-2018b"
  shell:
     """
        linenum0=`echo "($(wc -l < "{input[0]}") / 10) + 1" | bc`
        linenum1=`echo "($(wc -l < "{input[1]}") / 10) + 1" | bc`
        split -d {input[0]} {params.prefix[0]} -l $linenum0
        split -d {input[1]} {params.prefix[1]} -l $linenum1
        
     """

rule test:
  input:
    expand("splitbam/treat{index}_input{index}_allT_motif{{splitid}}",index=["1","2"])
  output:
    expand("site/treat{index}_input{index}_pvalue{{splitid}}",index=["1","2"])
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/calling_mRNA_tworep_clip.r {input[0]} {output[0]}
        Rscript code/calling_mRNA_tworep_clip.r {input[1]} {output[1]}
     """
