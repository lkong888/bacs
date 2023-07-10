#!/bin/bash

#SBATCH --job-name test     # Name for your job
#SBATCH --output mimseq.out    # Standard out goes to this file
#SBATCH --error mimseq.err     # Standard err goes to this file
#SBATCH -p short
#SBATCH --cpus-per-task=10

module load BEDTools/2.30.0-GCC-10.2.0
#zcat splitbam/treat1.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0' | awk '{OFS="\t"}{print $1,$2-2,$3+2,$4,$5,$6,$7,$8, $8/($7+$8)}' > splitbam/treat1.bed
#bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed splitbam/treat1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> splitbam/treat1.fasta.bed
#zcat splitbam/treat1.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0' | bedtools intersect -a - -b splitbam/treat1.fasta.bed -wa -wb -loj | awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > splitbam/treat1.sta.txt

#bedtools intersect -a <(cat splitbam/treat1.sta.txt |sort -k1,1 -k2,2n) -b <(zcat splitbam/input1.mRNA_pseusite.mpile.txt.gz| sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '($16+$17)>0'|awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$15,$16,$17,$17/($16+$17)}' > splitbam/treat1_input1_allT_motif.bed

#zcat splitbam/treat2.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0' | awk '{OFS="\t"}{print $1,$2-2,$3+2,$4,$5,$6,$7,$8, $8/($7+$8)}' > splitbam/treat2.bed
#bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed splitbam/treat2.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> splitbam/treat2.fasta.bed
#zcat splitbam/treat2.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0' | bedtools intersect -a - -b splitbam/treat2.fasta.bed -wa -wb -loj | awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > splitbam/treat2.sta.txt

#bedtools intersect -a <(cat splitbam/treat2.sta.txt | sort -k1,1 -k2,2n ) -b <(zcat splitbam/input2.mRNA_pseusite.mpile.txt.gz |  sort -k1,1 -k2,2n ) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '($16+$17)>0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$15,$16,$17,$17/($16+$17)}' > splitbam/treat2_input2_allT_motif.bed


####Rscript code/calling_mRNA1.r



#####treat1 VS input1
#zcat mergebam/treat1.mRNA_pseusite.mpile.txt.gz| awk '$7+$8>0' | awk '$8/($7+$8)>0'| awk '{OFS="\t"}{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}' > mergebam/treat1.bed
#bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed mergebam/treat1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> mergebam/treat1.fasta.bed
#zcat mergebam/treat1.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0'| bedtools intersect -a - -b mergebam/treat1.fasta.bed -wa -wb -loj | \
#        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}' > mergebam/treat1.sta.txt

#bedtools intersect -a <(cat mergebam/treat1.sta.txt |sort -k1,1 -k2,2n) -b <(zcat mergebam/input1.mRNA_pseusite.mpile.txt.gz| sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '$17+$18'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$9/$6,$10,$15,$16,$17,$18,$18/($17+$18),$19,$19/$16}' > mergebam/treat1_input1_allT_motif.bed

#colnames chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio

####treat2 vs input2
#zcat mergebam/treat2.mRNA_pseusite.mpile.txt.gz|awk '$7+$8>0' | awk '$8/($7+$8)>0' | awk '{OFS="\t"}{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}' > mergebam/treat2.bed
#bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed mergebam/treat2.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> mergebam/treat2.fasta.bed
#zcat mergebam/treat2.mRNA_pseusite.mpile.txt.gz |awk '$7+$8>0' | awk '$8/($7+$8)>0'| bedtools intersect -a - -b mergebam/treat2.fasta.bed -wa -wb -loj | \
#        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}' > mergebam/treat2.sta.txt

#bedtools intersect -a <(cat mergebam/treat2.sta.txt |sort -k1,1 -k2,2n) -b <(zcat mergebam/input2.mRNA_pseusite.mpile.txt.gz | sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '$17+$18'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$9/$6,$10,$15,$16,$17,$18,$18/($17+$18),$19,$19/$16}' > mergebam/treat2_input2_allT_motif.bed


module load R/4.0.3-foss-2020b
Rscript code/calling_mRNA.r
















