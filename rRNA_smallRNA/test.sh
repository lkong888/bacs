module load samtools/1.8-gcc5.4.0
module load Bowtie2/2.4.4-GCC-10.3.0

##########this script is to explore uniq alignment, best alignment and multiple alignment
####filter bam
###uniq alignment
module load SAMtools/0.1.19-GCC-10.3.0 ##total:23666714
samtools view test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep -v XS | cat <(samtools view -H test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam
samtools index test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam   ##4726512

##best alignment AS>XS
samtools view test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13>$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam
samtools index test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam ##10264304

###ambiguous alignment
samtools view test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13==$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam
samtools index test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam ##8675898

samtools view test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam | cut -f3 | sort | uniq -c > ambiguous.txt
samtools merge test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort_uniq_best.bam test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam 
samtools view test/input1.snoRNA.tRNA.clean.bowtie2.filter.sort_uniq_best.bam | cut -f3 | sort | uniq -c > determined_reads.txt

join -t $'\t' -1 2 -5 2 -e "0" -o auto -a1 -a2 <(cat determined_reads.txt | sed 's/^[ ]*//g' | sed 's/ /\t/g' |sort -k2) <(cat ambiguous.txt | sed 's/^[ ]*//g' | sed 's/ /\t/g' | sort -k2) > determined_ambiguous_reads.txt


####for treat1
samtools view test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep -v XS | cat <(samtools view -H test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam
samtools index test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam   ##4726512

##best alignment AS>XS
samtools view test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13>$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam
samtools index test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam ##10264304

###ambiguous alignment
samtools view test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13==$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam
samtools index test/treat1.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam ##8675898



samtools view test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep -v XS | cat <(samtools view -H test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam
samtools index test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.uniq.bam   ##4726512

##best alignment AS>XS
samtools view test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13>$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam
samtools index test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.best.bam ##10264304

###ambiguous alignment
samtools view test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam | grep XS | sed 's/AS:i:/AS:i\t/g' | sed 's/XS:i:/XS:i\t/g' | awk '$13==$15' |\
 sed 's/AS:i\t/AS:i:/g' | sed 's/XS:i\t/XS:i:/g' | cat <(samtools view -H test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.bam) - | samtools view -bS - > test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam
samtools index test/treat2.snoRNA.tRNA.clean.bowtie2.filter.sort.ambiguous.bam ##8675898





##################for using only uniq and best alignment(AS>XS)
cat test/treat1.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat1.uniq.best.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat1.uniq.best.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat1.uniq.best.fasta.bed
cat test/treat1.uniq.best.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat1.uniq.best.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat1.uniq.best.sta.txt

cat test/input1.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input1.uniq.best.bed
bedtools intersect -a test/treat1.uniq.best.sta.txt -b test/input1.uniq.best.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat1_input1.uniq.best_motif.bed

#treat2 vs input2
cat test/treat2.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat2.uniq.best.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat2.uniq.best.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat2.uniq.best.fasta.bed
cat test/treat2.uniq.best.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat2.uniq.best.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat2.uniq.best.sta.txt

cat test/input2.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.uniq.best.bed
bedtools intersect -a test/treat2.uniq.best.sta.txt -b test/input2.uniq.best.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat2_input2.uniq.best_motif.bed

#treat3 vs input2
cat test/treat3.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat3.uniq.best.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat3.uniq.best.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat3.uniq.best.fasta.bed
cat test/treat3.uniq.best.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat3.uniq.best.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat3.uniq.best.sta.txt

cat test/input2.snoRNA.tRNA.filter.uniq.best_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.uniq.best.bed
bedtools intersect -a test/treat3.uniq.best.sta.txt -b test/input2.uniq.best.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat3_input2.uniq.best_motif.bed





#########################find T to A /G higher than T ro C in snoRNA and tRNA
cat treat2.snoRNA.tRNA.filter.uniq.best.mpile.txt | awk '$3=="T" || $3=="t"' | awk '$4>20' |awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5+$14, $6+$15, $7+$16, ($5+$14)/($4), ($6+$15)/$4, ($7+$16)/$4}}' | awk '$9> 0.02' | awk '($8+$10)>$9' > treat2_snoRNA_tRNA_exclude.txt
cat treat3.snoRNA.tRNA.filter.uniq.best.mpile.txt | awk '$3=="T" || $3=="t"' | awk '$4>20' |awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5+$14, $6+$15, $7+$16, ($5+$14)/($4), ($6+$15)/$4, ($7+$16)/$4}}' | awk '$9> 0.02' | awk '($8+$10)>$9' > treat3_snoRNA_tRNA_exclude.txt



cat treat2.snoRNA.tRNA.filter.uniq.best.mpile.txt | awk '$3=="T" || $3=="t"' | awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20, $5+$14+$7+$16}}' > treat2.snoRNA.tRNA.filter.uniq.best.mpile_T.txt
cat treat3.snoRNA.tRNA.filter.uniq.best.mpile.txt | awk '$3=="T" || $3=="t"' | awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $8+$17, $6+$15, $11+$20, $5+$14+$7+$16}}' > treat3.snoRNA.tRNA.filter.uniq.best.mpile_T.txt

 

###compare the sites using all reads and using only uniq and best alignments
cat snoRNA_tRNA.uniq.best_called_sites.txt | sed 's/ /\t/g' | wc -l ##1089
cat snoRNA_tRNA_called_sites.txt  | sed 's/ /\t/g' | wc -l #1170

bedtools intersect -a <(cat snoRNA_tRNA.uniq.best_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.uniq.best_called_sites_only.txt ##53
bedtools intersect -a <(cat snoRNA_tRNA_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.uniq.best_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA_called_sites_only.txt ##134
bedtools intersect -a <(cat snoRNA_tRNA_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.uniq.best_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) > snoRNA_tRNA_called_sites_overlap.txt ##1035



##################for using only mapq 10
cat test/treat1.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat1.mapq10.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat1.mapq10.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat1.mapq10.fasta.bed
cat test/treat1.mapq10.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat1.mapq10.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat1.mapq10.sta.txt

cat test/input1.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input1.mapq10.bed
bedtools intersect -a test/treat1.mapq10.sta.txt -b test/input1.mapq10.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat1_input1.mapq10_motif.bed

#treat2 vs input2
cat test/treat2.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat2.mapq10.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat2.mapq10.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat2.mapq10.fasta.bed
cat test/treat2.mapq10.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat2.mapq10.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat2.mapq10.sta.txt

cat test/input2.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.mapq10.bed
bedtools intersect -a test/treat2.mapq10.sta.txt -b test/input2.mapq10.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat2_input2.mapq10_motif.bed

#treat3 vs input2
cat test/treat3.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat3.mapq10.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat3.mapq10.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat3.mapq10.fasta.bed
cat test/treat3.mapq10.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat3.mapq10.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat3.mapq10.sta.txt

cat test/input2.snoRNA.tRNA.filter.mapq10_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.mapq10.bed
bedtools intersect -a test/treat3.mapq10.sta.txt -b test/input2.mapq10.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat3_input2.mapq10_motif.bed







##################for using only mapq 1
cat test/treat1.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat1.mapq1.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat1.mapq1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat1.mapq1.fasta.bed
cat test/treat1.mapq1.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat1.mapq1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat1.mapq1.sta.txt

cat test/input1.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input1.mapq1.bed
bedtools intersect -a test/treat1.mapq1.sta.txt -b test/input1.mapq1.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat1_input1.mapq1_motif.bed

#treat2 vs input2
cat test/treat2.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat2.mapq1.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat2.mapq1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat2.mapq1.fasta.bed
cat test/treat2.mapq1.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat2.mapq1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat2.mapq1.sta.txt

cat test/input2.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.mapq1.bed
bedtools intersect -a test/treat2.mapq1.sta.txt -b test/input2.mapq1.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat2_input2.mapq1_motif.bed

#treat3 vs input2
cat test/treat3.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat3.mapq1.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat3.mapq1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat3.mapq1.fasta.bed
cat test/treat3.mapq1.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat3.mapq1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat3.mapq1.sta.txt

cat test/input2.snoRNA.tRNA.filter.mapq1_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.mapq1.bed
bedtools intersect -a test/treat3.mapq1.sta.txt -b test/input2.mapq1.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat3_input2.mapq1_motif.bed





####compare using all reads and using only mapq1
bedtools intersect -a <(cat snoRNA_tRNA.mapq1_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.mapq1_only_called_sites.txt ##49
bedtools intersect -a <(cat snoRNA_tRNA_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.mapq1_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA_called_sites_only_compare_mapq1.txt ##25


bedtools intersect -a <(cat snoRNA_tRNA.mapq10_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.mapq10_only_called_sites.txt ##49
bedtools intersect -a <(cat snoRNA_tRNA_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.mapq10_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA_called_sites_only_compare_mapq10.txt ##25



####compare using AS>XS reads and using only mapq1 or mapq 10
bedtools intersect -a <(cat snoRNA_tRNA.mapq1_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.uniq.best_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.mapq1_only_called_sites_against_best.txt ##49
bedtools intersect -a <(cat snoRNA_tRNA.uniq.best_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.mapq1_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.uniq.best_called_sites_only_compare_mapq1.txt ##25


bedtools intersect -a <(cat snoRNA_tRNA.mapq10_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.uniq.best_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.mapq10_only_called_sites_against_best.txt ##44
bedtools intersect -a <(cat snoRNA_tRNA.uniq.best_called_sites.txt| sed 's/ /\t/g' | grep -v ^chr) -b <(cat snoRNA_tRNA.mapq10_called_sites.txt | sed 's/ /\t/g' | grep -v ^chr) -v > snoRNA_tRNA.uniq.best_called_sites_only_compare_mapq10.txt ##182







#######################################################upstream and downstream 10 bases
cat align/treat1.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat1.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat1.fasta.bed
cat test/treat1.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat1.sta.txt

cat align/input1.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input1.bed
bedtools intersect -a test/treat1.sta.txt -b test/input1.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat1_input1_motif.bed

#treat2 vs input2
cat align/treat2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat2.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat2.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat2.fasta.bed
cat test/treat2.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat2.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat2.sta.txt

cat align/input2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.bed
bedtools intersect -a test/treat2.sta.txt -b test/input2.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat2_input2_motif.bed

#treat3 vs input2
cat align/treat3.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>6' | awk '{OFS="\t"}{print $1, $2-6,$2+5, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > test/treat3.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/treat3.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/treat3.fasta.bed
cat test/treat3.bed | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b test/treat3.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > test/treat3.sta.txt

cat align/input2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > test/input2.bed
bedtools intersect -a test/treat3.sta.txt -b test/input2.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > test/treat3_input2_motif.bed



########get up and down stream 10 bp

/users/ludwig/ebu571/ebu571/20May2023_final/site/snoRNA_tRNA_called_sites_fdrpassed.txt

cat site/snoRNA_tRNA_called_sites_fdrpassed.txt| grep -v ^chr | awk '$2>6' | awk '{OFS="\t"}{print $1, $2-5,$3+5,$4,$5,$6,$7}' |sed 's/-/_/g' > test/snoRNA_tRNA_called_sites_fdrpassed.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/snoRNA_tRNA_called_sites_fdrpassed.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+5, $3-5, $4}'> test/snoRNA_tRNA_called_sites_fdrpassed.fasta.bed
cat site/snoRNA_tRNA_called_sites_fdrpassed.txt| grep -v ^chr | bedtools intersect -a - -b test/snoRNA_tRNA_called_sites_fdrpassed.fasta.bed -wa -wb -loj > test/snoRNA_tRNA_called_sites_fdrpassed.motif.bed



#######find motif
module load BEDTools/2.30.0-GCC-10.2.0
cat site/snoRNA_tRNA_called_sites_fdrpassed.txt| grep -v ^chr | grep -v tRNA | awk '$2>6' | awk '{OFS="\t"}{print $1, $2-5,$3+5,$4,$5,$6,$7}' |sed 's/-/_/g' > test/snoRNA_called_sites_fdrpassed.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed test/snoRNA_called_sites_fdrpassed.bed > test/snoRNA_called_sites_fdrpassed.fasta
/well/ludwig/users/ebu571/conda/skylake/envs/meme/bin/meme test/snoRNA_called_sites_fdrpassed.fasta -rna 


#######conda install weblogo
module load Anaconda3/2021.05
eval "$(conda shell.bash hook)"
conda create -n weblogo
conda activate weblogo
conda install -c conda-forge weblogo


/well/ludwig/users/ebu571/conda/skylake/envs/meme/bin/meme test/snoRNA_called_sites_fdrpassed.fasta -rna 