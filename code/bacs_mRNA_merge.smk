"""
Workflow to compute taps conversion(notrim) #####merge bam file of both 26May2023 and 23June2023
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 8 --snakefile code/bacs_mRNA_merge.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs1 --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["merge"]
print(SAMPLES)

rule all:
    input:
        #expand("mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam", sample = SAMPLES),
        #expand("mergebam/{sample}.mRNA_first.mpile.txt.gz", sample = SAMPLES),
       # expand("mergebam/{sample}.mapping_report.txt", sample = SAMPLES),
        #expand("mergebam/{sample}.mRNA_pseusite.mpile.txt.gz", sample = SAMPLES),
        #expand("rpkm/{sample}-chrM.tpm", sample = SAMPLES),
        expand("mergebam/treat{index}_input{index}_allT_motif{splitid}", index=['1','2'], splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("site/treat{index}_input{index}_pvalue{splitid}", index=['1','2'], splitid=['00','01','02','03','04','05','06','07','08','09']),
        expand("site/treat{index}_input{index}_site.txt", index=['1','2'])
#######################################################################merge mRNA#############################################################
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


rule strand_clasify: 
    input:
        "mergebam/{sample}.mrna.filter_dedup.sort.bam"
    output:
        first="mergebam/{sample}.mRNA_first.bam",
        second="mergebam/{sample}.mRNA_second.bam",
    params:
        mq=" -q 30"
    envmodules: 
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -h {input} |awk '$2==0 ||$1~/^@/' |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - >{output.first}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -h {input} |awk '$2==16 ||$1~/^@/' |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - >{output.second}
        """

rule count:
    input:
        first="mergebam/{sample}.mRNA_first.bam",
        second="mergebam/{sample}.mRNA_second.bam",
    output:
        first="mergebam/{sample}.mRNA_first.mpile.txt.gz",
        second="mergebam/{sample}.mRNA_second.mpile.txt.gz"
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
        first="mergebam/{sample}.mRNA_first.mpile.txt.gz",
        second="mergebam/{sample}.mRNA_second.mpile.txt.gz"
    output:
        "mergebam/{sample}.mRNA_pseusite.mpile.txt.gz"
    shell:
        """
        zcat {input.first}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>0'| awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "-", $3, $4, $5+$14, $7+$16,$11+$20}}' | gzip > {output}
        zcat {input.second}  | grep ^chr | sed 's/,/\\t/g' | awk '$4>0' | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2-1, $2, "+", $3, $4, $8+$17, $6+$15, $11+$20}}' | gzip >> {output}

        """


###chr     pos     ref_base        depth,A,C,G,T,N,Skip,Gap,Insert,Delete,a,c,g,t,n,skip,gap,insert,delete
###chr     pos     ref_base      strand  depth A G gap
###chr     pos     ref_base      strand  depth T C gap


###cat H41_2023Feb14_S40.mRNA_pseusite.mpile.txt | awk '$7+$8>5' | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $8/($7+$8)}' | awk '$9>=0.05' | wc -l
##22361

#cat H43_2023Feb14_S42.mRNA_pseusite.mpile.txt | awk '$7+$8>0' | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $8/($7+$8)}' | awk '$9>=0.01' > H43_ctrl.bed
##2467

##cat H41_2023Feb14_S40.mRNA_pseusite.mpile.txt | awk '$7+$8>5' | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $8/($7+$8)}' | awk '$9>=0.05' | bedtools intersect -a - -b H43_ctrl.bed -v -wa > H41_mutated.bed
#cat H42_2023Feb14_S41.mRNA_pseusite.mpile.txt | awk '$7+$8>5' | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $8/($7+$8)}' | awk '$9>=0.05' | bedtools intersect -a - -b H43_ctrl.bed -v -wa > H42_mutated.bed
#cat H40_2023Feb14_S39.mRNA_pseusite.mpile.txt | awk '$7+$8>5' | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $8/($7+$8)}' | awk '$9>=0.05' | bedtools intersect -a - -b H43_ctrl.bed -v -wa > H40_mutated.bed

rule mapping:
    input:
        mergebam="mergebam/{sample}.mRNAAligned.sortedByCoord.out.bam",
        filt_dedup="mergebam/{sample}.mrna.filter_dedup.sort.bam"
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


##############call site
rule splitbed:
  input:
    expand("mergebam/treat{index}_input{index}_allT_motif.bed", index=['1','2'])
  output:
    expand("mergebam/treat{index}_input{index}_allT_motif{splitid}", index=['1','2'], splitid=['00','01','02','03','04','05','06','07','08','09'])
  params:
    prefix=expand("mergebam/treat{index}_input{index}_allT_motif",index=['1','2'])
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
    expand("mergebam/treat{index}_input{index}_allT_motif{{splitid}}", index=['1','2'])
  output:
    expand("site/treat{index}_input{index}_pvalue{{splitid}}", index=['1','2'])
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        Rscript code/calling_mRNA.r {input[0]} {output[0]}
        Rscript code/calling_mRNA.r {input[1]} {output[1]}
     """

rule padjust:
  input:
    expand("site/treat{{index}}_input{{index}}_pvalue{splitid}",splitid=['00','01','02','03','04','05','06','07','08','09'])
  output:
    tmp=expand("site/treat{{index}}_input{{index}}_site.tmp.txt"),
    site=expand("site/treat{{index}}_input{{index}}_site.txt")
  threads: 1
  envmodules:
    "R/4.0.3-foss-2020b"
  shell:
     """
        cat {input} | grep -v ^motif | awk '{{OFS="\\t"}}{{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}}' > {output.tmp}
        Rscript code/fdr.r {output.tmp} {output.site}
        rm {output.tmp}
     """