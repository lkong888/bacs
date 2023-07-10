"""
Workflow to compute taps conversion(notrim)
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 14 --snakefile code/analysis_for_haiqi.smk --cluster "sbatch -p short -o bacs.log -e bacs.err --job-name=bacs --constraint=hga --cpus-per-task 3 "
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["ribo_rna"]
nn_fa = config["nn_fa"]
nn_fasta = config["nn_fasta"]

print(SAMPLES)

rule all:
    input:
        expand("fastq/{sample}.clean.merge.fq.gz", sample = SAMPLES),
        expand("align/{sample}.nnunn.clean.filter.sort_pseusite.mpile.txt", sample = SAMPLES),
        expand("align/{sample}.rRNA.filter.dedup.bam", sample = SAMPLES),
        expand("align/{sample}.rRNA.filter.dedup_pseusite.mpile.txt", sample = SAMPLES),
        expand("align/{sample}.clean.nnunn.sort_NNUNN1.contex.table", sample = SAMPLES),
        expand("align/{sample}.rRNA.filter.dedup.mpile.2kb_mutation_sta_plot.txt", sample = SAMPLES),
        expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam", sample = SAMPLES),
        expand("align/{sample}.mapping_report.txt", sample = SAMPLES),
        expand("align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam", sample = SAMPLES)



##############################################################################################################################################################################
############################reads cleaning up ########################

rule rm_adapter:
    input:
        expand("fastq/{{sample}}_R{readDirection}.fastq.gz",readDirection=['1','2'])
    output:
        cut1=expand("fastq/{{sample}}.R{readDirection}.trim1.fastq.gz", readDirection=['1','2']),
        cut2=expand("fastq/{{sample}}.R{readDirection}.trim.fq.gz", readDirection=['1','2'])
    params:
        prefix="fastq/{sample}.cutadapt"
    threads: 3
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/cutadapt/bin/cutadapt -j 4 -e 0.15 -n 1 -O 3 -q 20 -m 18 \
         -g GACGCTCTTCCGATCT -A AGATCGGAAGAGCGTC -o {output.cut1[0]} -p {output.cut1[1]} \
         {input} > {params.prefix}_round1.log

         /well/ludwig/users/ebu571/conda/skylake/envs/cutadapt/bin/cutadapt -j 4 -e 0.15 -n 1 -O 9 -q 20 -m 18 --nextseq-trim=20 \
         -a NNNNNNAGATCGGAAGAGCACA -G TGTGCTCTTCCGATCT -o {output.cut2[0]} -p {output.cut2[1]} \
        {output.cut1[0]} {output.cut1[1]} > {params.prefix}_round2.log
        """
rule rm_UMI:
    input:
        expand("fastq/{{sample}}.R{readDirection}.trim.fq.gz", readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}.R{readDirection}.clean.fq.gz", readDirection=['1','2'])
    threads: 3
    params:
        prefix="fastq/{sample}.deUMI"
    envmodules: "UMI-tools/1.0.1-foss-2019b-Python-3.7.4"
    shell:
        """
        umi_tools extract --random-seed=123 --extract-method=regex --bc-pattern=".*" --bc-pattern2="^(?P<umi_1>.{{6}}).*" -I {input[0]} --read2-in={input[1]} \
	    --stdout={output[0]} --read2-out={output[1]} -L {params.prefix}.log
        """
###############################################join paired reads ##############################################################################
rule merge:
    input:
        expand("fastq/{{sample}}.R{readDirection}.clean.fq.gz", readDirection=['1','2'])
    output:
        m="fastq/{sample}.clean.merge.fq.gz",
        u1="merge_reads/{sample}.clean.pass.unmerge_r1.fq.gz",
        u2="merge_reads/{sample}.clean.pass.unmerge_r2.fq.gz",
        f1="merge_reads/{sample}.clean.r1passr2fail.fq.gz",
        f2="merge_reads/{sample}.clean.r1failr2pass.fq.gz",
    log:
        "logs/{sample}.trim.log"
    params:
        "merge_reads/{sample}.merge"
    threads: 3
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/fastp/bin/fastp --thread {threads} -i {input[0]} \
            --disable_adapter_trimming --merge --correction --overlap_len_require 10 --overlap_diff_percent_limit 20 \
            -i {input[0]} -I {input[1]} --merged_out {output.m} --out1 {output.u1} --out2 {output.u2} --unpaired1 {output.f1} --unpaired2 {output.f2} -h {params}.html -j {params}.json
        """
############################Align for NNUNN spikein########################
rule align_nnunn: 
    input:
        "fastq/{sample}.clean.merge.fq.gz",
    output:
        bam="align/{sample}.clean.nnunn.bam",
        unalign="fastq/{sample}.nnunn.clean_unalign.fq.gz"
    log:
        "logs/{sample}.clean.nnunn.log"
    params:
        tmp="fastq/{sample}.nnunn.clean_unalign.fq.gz",
        ref=nn_fa
    envmodules: "Bowtie2/2.4.4-GCC-10.3.0"
    threads: 2
    shell:
        """
        (bowtie2 -p 2 --no-unal --local --score-min G,20,1 --n-ceil L,0,0.25 -N 1 --mp 3 -L 10 --np 0 --un-gz {params.tmp} -x {params.ref} -U {input} | \
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -o {output.bam} )  1>{log} 2>&1

        """

rule nnunn_filter:
    input:
        align="align/{sample}.clean.nnunn.bam",
    output:
        filter_bam="align/{sample}.clean.nnunn.filter.bam",
        disgard_bam="align/{sample}.clean.nnunn.discard.bam",
        sortbam="align/{sample}.clean.nnunn.sort.bam",
        fqgz="fastq/{sample}.nnunn.nnunn_discard.fq.gz"
    params:
        mq=" -q 10",
        unalifq="fastq/{sample}.nnunn.nnunn_discard.fq",
        tmp="{sample}.tmp",
    envmodules:"BEDTools/2.30.0-GCC-10.2.0"
    shell:
        """
        join -t $'\\t' -1 1 -2 1 <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {input.align} | awk '{{OFS="\\t"}}{{print $0, length($10)}}' | sort -k 1) <(bedtools bamtobed -i {input.align} | \
            awk '($3-$2)>0' | sed 's/\/.//g' | awk '{{OFS="\\t"}}{{print $4, $3-$2}}' | sort -k 1) | \
            awk '$NF/$(NF-1)>=0.9' | awk '{{OFS="\\t"}}{{for(i=0;++i<=NF-3;) printf $i"\\t"; print $(NF-2)}}' |\
            cat <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view  {input.align} -H) - |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view - -bS >> {output.filter_bam}

        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {params.mq} {output.filter_bam} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=10' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}

        join -t $'\\t' -1 1 -2 1 <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {input.align} | awk '{{OFS="\\t"}}{{print $0, length($10)}}' | sort -k 1) <(bedtools bamtobed -i {input.align} | \
            awk '($3-$2)>0' | sed 's/\/.//g' | awk '{{OFS="\\t"}}{{print $4, $3-$2}}' | sort -k 1) | \
            awk '$NF/$(NF-1)<0.9' | awk '{{OFS="\\t"}}{{for(i=0;++i<=NF-3;) printf $i"\\t"; print $(NF-2)}}' |\
            cat <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view  {input.align} -H) - |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view - -bS >> {output.disgard_bam}

        bedtools bamtofastq -i {output.disgard_bam} \
                      -fq {params.unalifq} 
        gzip {params.unalifq}
        """
##only keep reads that have (mapped length)/(read length)>=0.9
##for reads that mapped but with (mapped length)/(read length)<0.9, conver to fastq and merge into unaligned reads

rule nnunn_count:
    input:
        "align/{sample}.clean.nnunn.sort.bam"
    output:
        mpile="align/{sample}.nnunn.clean.filter.sort.mpile.txt",
        sta="align/{sample}.nnunn.clean.filter.sort_pseusite.mpile.txt"
    params:
        ref=nn_fasta
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {input} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15}}' > {output.sta}
        """

rule nnunn_context:
    input:
        "align/{sample}.clean.nnunn.sort.bam"
    output:
        "align/{sample}.clean.nnunn.sort.bam.bai",
        "align/{sample}.clean.nnunn.sort_NNUNN1.contex.table"
    envmodules: "R/4.0.3-foss-2020b "
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {input}
        Rscript code/nnunn_sta.r {input}
        """

##############################################################################################################################################################################
############################Align for rRNA spikein########################       
rule align_rRNA: 
    input:
        unalign="fastq/{sample}.nnunn.clean_unalign.fq.gz",
        discard="fastq/{sample}.nnunn.nnunn_discard.fq.gz"
    output:
        unmap="fastq/{sample}.nnunn.clean_unmap.fq.gz",
        bam="align/{sample}.rRNA.clean.bowtie2.bam",
        unalign="fastq/{sample}.rRNA.clean_unalign.fq.gz"
    log:
        "logs/{sample}.rRNA.clean.bowtie2.log"
    params:
        tmp="fastq/{sample}.rRNA.clean_unalign.fq.gz",
        ref="resource/spikein_rRNA"
    envmodules: "Bowtie2/2.4.4-GCC-10.3.0"
    threads: 2
    shell:
        """
        zcat {input.unalign} {input.discard} > {output.unmap}
        (bowtie2 -p 2 --no-unal --local -L 16 -N 1 --mp 4 \
        --un-gz {params.tmp} -x {params.ref} -U {output.unmap} | \
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -o {output.bam} ) 1>{log} 2>&1
        
        """
rule filter_split: 
    input:
        "align/{sample}.rRNA.clean.bowtie2.bam"
    output:
        sortbam="align/{sample}.rRNA.clean.bowtie2.filter.sort.bam",
        longbam="align/{sample}.rRNA.clean.bowtie2.filter_long.sort.bam",
        shortbam="align/{sample}.rRNA.clean.bowtie2.filter_short.sort.bam"
    params:
        mq=" -q 10",
        tmp="{sample}.tmp",
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {params.mq} {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=10' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}
        
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {output.sortbam} NR_003286.4_RNA18SN5 NR_003287.4_RNA28SN5 MT_RNR1_ENST00000389680.2_ENSE00001544499 MT_RNR2_ENST00000387347.2_ENSE00001544497 >{output.longbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.longbam}

        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {output.sortbam} | awk '$3!="NR_003286.4_RNA18SN5" && $3!="NR_003287.4_RNA28SN5" && $3!="MT_RNR1_ENST00000389680.2_ENSE00001544499" && $3!="MT_RNR2_ENST00000387347.2_ENSE00001544497"' | \
        cat <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {output.sortbam} -H) - |/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view - -bS > {output.shortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.shortbam}
        """

rule rrna_dedup: 
    input:
        "align/{sample}.rRNA.clean.bowtie2.filter_long.sort.bam"
    output:
        dedup="align/{sample}.rRNA.clean.bowtie2.filter_long.sort_dedup.bam"
    params:
        mq=" -q 10",
        tmp="{sample}.tmp",
        prefix="logs/{sample}.rrna.dedup"
    envmodules: "UMI-tools/1.0.1-foss-2019b-Python-3.7.4"
    threads: 2
    shell:
        """
        umi_tools dedup --random-seed=123 --method=unique -I {input} -S {output.dedup} --output-stats {params.prefix}.stats -L {params.prefix}.logs
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.dedup}
        """

rule merge_count: 
    input:
        longbam="align/{sample}.rRNA.clean.bowtie2.filter_long.sort_dedup.bam",
        shortbam="align/{sample}.rRNA.clean.bowtie2.filter_short.sort.bam"
    output:
        merge="align/{sample}.rRNA.filter.dedup.bam",
        mpile="align/{sample}.rRNA.filter.dedup.mpile.txt",
        sta="align/{sample}.rRNA.filter.dedup_pseusite.mpile.txt"
    params:
        ref="resource/spikein_rRNA.fa"
    envmodules: "CMake/3.23.1-GCCcore-11.3.0", "R/4.2.1-foss-2022a "
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.merge} {input.longbam} {input.shortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools mpileup -d 0 -Q 10 --reverse-del \
            -f {params.ref} {output.merge} | /users/ludwig/ebu571/ebu571/tools/cpup/cpup | sed 's/,/\\t/g' > {output.mpile}
        cat {output.mpile} | grep ^spike_in | awk '$3=="A" || $3=="a"' | awk '{{OFS="\\t"}}{{print $1, $2,  $3, $4, $5+$14, $7+$16, $11+$20}}' > {output.sta}
        cat {output.mpile} | grep -v ^spike_in | grep -v ^nnunn| awk '$3=="T" || $3=="t"' | awk '{{OFS="\\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20}}' >> {output.sta}
        
        """

rule unmodified_2kb:
    input:
        mpile="align/{sample}.rRNA.filter.dedup.mpile.txt"
    output:
        "align/{sample}.rRNA.filter.dedup.mpile.2kb_mutation_sta_plot.txt"
    envmodules: "R/4.2.1-foss-2022a "
    threads: 2
    shell:
        """
            Rscript code/spikein_sta_2kb.r {input} 
        """
##############################################################################################################################################################################
############################Align for smallRNA spikein########################     
rule align_snorna: 
    input:
        "fastq/{sample}.rRNA.clean_unalign.fq.gz"
    output:
        bam="align/{sample}.snoRNA.clean.bowtie2.bam",
        unalign="fastq/{sample}.snoRNA.clean_unalign.fq.gz"
    log:
        "logs/{sample}.snoRNA.clean.bowtie2.log"
    params:
        tmp="fastq/{sample}.snoRNA.clean_unalign.fq.gz",
        ref="resource/snoRNA_reference"
    envmodules: "Bowtie2/2.4.4-GCC-10.3.0"
    threads: 2
    shell:
        """
        (bowtie2 -p 2 --no-unal --local -L 16 -N 1 --mp 4\
        --un-gz {params.tmp} -x {params.ref} -U {input} | \
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -o {output.bam} ) 1>{log} 2>&1
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.bam}
        """

##############################################################################################################################################################################
############################Align for tRNA spikein########################     
rule align_trna: 
    input:
        "fastq/{sample}.snoRNA.clean_unalign.fq.gz"
    output:
        bam="align/{sample}.tRNA.clean.bowtie2.bam",
        unalign="fastq/{sample}.tRNA.clean_unalign.fq.gz"
    log:
        "logs/{sample}.tRNA.clean.bowtie2.log"
    params:
        tmp="fastq/{sample}.tRNA.clean_unalign.fq.gz",
        ref="resource/tRNA_high_confident"
    envmodules: "Bowtie2/2.4.4-GCC-10.3.0"
    threads: 2
    shell:
        """
        (bowtie2 -p 2 --no-unal --local -L 16 -N 1 --mp 4\
        --un-gz {params.tmp} -x {params.ref} -U {input}| \
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -o {output.bam} ) 1>{log} 2>&1
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.bam}
        """


rule merge_sno_trna: 
    input:
        sno="align/{sample}.snoRNA.clean.bowtie2.bam",
        trna="align/{sample}.tRNA.clean.bowtie2.bam"

    output:
        "align/{sample}.snoRNA.tRNA.clean.bowtie2.bam"
        
    envmodules: "samtools/1.8-gcc5.4.0"
    threads: 2
    shell:
        """
        samtools merge {output} {input.sno} {input.trna}
        """
rule filter: 
    input:
        "align/{sample}.snoRNA.tRNA.clean.bowtie2.bam"
    output:
        sortbam="align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam",
    params:
        tmp="{sample}.tmp",
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {input}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS {input} | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools sort -@ 8 -O BAM -T {params.tmp} --input-fmt-option 'filter=[NM]<=10' >{output.sortbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools index {output.sortbam}
        """
rule best: 
    input:
        "align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam"
    output:
        uniqbam="align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.bam",
        bestbam="align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.best.bam",
        bam="align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.uniq.best.bam"
    params:
        tmp="{sample}.tmp",
    threads: 2
    shell:
        """
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {input} | grep -v XS | cat <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -H {input}) - | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - > {output.uniqbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view {input} | grep XS | sed 's/AS:i:/AS:i\\t/g' | sed 's/XS:i:/XS:i\\t/g' | awk '$13>$15' | \
                sed 's/AS:i\\t/AS:i:/g' | sed 's/XS:i\\t/XS:i:/g' | cat <(/well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -H {input}) - | /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools view -bS - > {output.bestbam}
        /well/ludwig/users/ebu571/conda/skylake/envs/samtools/bin/samtools merge {output.bam} {output.uniqbam} {output.bestbam}

        """

rule mapping:
    input:
        fq1="fastq/{sample}_R1.fastq.gz",
        fq2="fastq/{sample}_R2.fastq.gz",
        clean1="fastq/{sample}.R1.clean.fq.gz",
        clean2="fastq/{sample}.R2.clean.fq.gz",
        mergefq="fastq/{sample}.clean.merge.fq.gz",
        nnunn="align/{sample}.clean.nnunn.filter.bam",
        nnunn_q10="align/{sample}.clean.nnunn.sort.bam",
        rrna_input="fastq/{sample}.nnunn.clean_unmap.fq.gz",
        rrna_mapped="align/{sample}.rRNA.clean.bowtie2.bam",
        rrna_filtered="align/{sample}.rRNA.filter.dedup.bam",
        rrna_unalign="fastq/{sample}.rRNA.clean_unalign.fq.gz",
        snorna_map="align/{sample}.snoRNA.clean.bowtie2.bam",
        snorna_unalign="fastq/{sample}.snoRNA.clean_unalign.fq.gz",
        trna_map="align/{sample}.tRNA.clean.bowtie2.bam",
        trna_unalign="fastq/{sample}.tRNA.clean_unalign.fq.gz",
        snorna_trna_merge="align/{sample}.snoRNA.tRNA.clean.bowtie2.bam",
        snorna_trna_filtered="align/{sample}.snoRNA.tRNA.clean.bowtie2.filter.sort.bam",
    output:
        "align/{sample}.mapping_report.txt", 
    envmodules: "samtools/1.8-gcc5.4.0"
    shell:
        """
        echo -e "smp\\traw_reads_r1\\traw_reads_r2\\tclean_reads_r1\\tclean_reads_r2\\tmergefq\\tnnunn\\tnnunn_q10\\trrna_input\\trrna_mapped\\trrna_filtered\\trrna_unalign\\tsnorna_map\\tsnorna_unalign\\ttrna_map\\ttrna_unalign\\tsnorna_trna_merge\\tsnorna_trna_filtered" > {output}
        raw_reads_r1=`echo $(zcat {input.fq1}|wc -l)/4|bc`
        raw_reads_r2=`echo $(zcat {input.fq2}|wc -l)/4|bc`
        clean_reads_r1=`echo $(zcat {input.clean1}|wc -l)/4|bc`
        clean_reads_r2=`echo $(zcat {input.clean2}|wc -l)/4|bc`
        mergefq=`echo $(zcat {input.mergefq}|wc -l)/4|bc`
        nnunn=`samtools view {input.nnunn} | cut -f 1 | sort | uniq | wc -l`
        nnunn_q10=`samtools view {input.nnunn_q10} | cut -f 1 | sort | uniq | wc -l`
        rrna_input=`echo $(cat {input.rrna_input}|wc -l)/4|bc`
        rrna_mapped=`samtools view {input.rrna_mapped} | cut -f 1 | sort | uniq | wc -l`
        rrna_filtered=`samtools view {input.rrna_filtered} | cut -f 1 | sort | uniq | wc -l`
        rrna_unalign=`echo $(zcat {input.rrna_unalign}|wc -l)/4|bc`

        snorna_map=`samtools view {input.snorna_map} | cut -f 1 | sort | uniq | wc -l`
        snorna_unalign=`echo $(zcat {input.snorna_unalign}|wc -l)/4|bc`

        trna_map=`samtools view {input.trna_map} | cut -f 1 | sort | uniq | wc -l`
        trna_unalign=`echo $(zcat {input.trna_unalign}|wc -l)/4|bc`

        snorna_trna_merge=`samtools view {input.snorna_trna_merge} | cut -f 1 | sort | uniq | wc -l`
        snorna_trna_filtered=`samtools view {input.snorna_trna_filtered} | cut -f 1 | sort | uniq | wc -l`


        echo -e "{output}\\t${{raw_reads_r1}}\\t${{raw_reads_r2}}\\t${{clean_reads_r1}}\\t${{clean_reads_r2}}\\t${{mergefq}}\\t${{nnunn}}\\t${{nnunn_q10}}\\t${{rrna_input}}\\t${{rrna_mapped}}\\t${{rrna_filtered}}\\t${{rrna_unalign}}\\t${{snorna_map}}\\t${{snorna_unalign}}\\t${{trna_map}}\\t${{trna_unalign}}\\t${{snorna_trna_merge}}\\t${{snorna_trna_filtered}}" >> {output}
        
        """
