module load samtools/1.8-gcc5.4.0
module load Bowtie2/2.4.4-GCC-10.3.0

bowtie2-build --threads 2 nnunn.fa nnunn
bowtie2-build --threads 2 spikein_rRNA.fa spikein_rRNA
bowtie2-build --threads 2 tRNA_snorna_cdhit.fa tRNA_snorna_cdhit


###check conversion and false positives

echo -e "smp\tspikein29_ratio\tspikein29_depth\tspikein43_ratio\tspikein43_depth\tfp_ratio\tfp_depth" >all_conversion_fp.txt
for f in `ls *rRNA.filter.dedup_pseusite.mpile.txt` #`ls H*rRNA.snoRNA.tRNA_merge.filter.sort_pseusite.mpile.txt` 
do
    spikein29_ratio=`cat $f | awk '$1=="spike_in" && ($2=="29")' | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print C/(C+T)}'`
    spikein29_depth=`cat $f | awk '$1=="spike_in" && ($2=="29")' | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print (C+T)}'`
    spikein43_ratio=`cat $f | awk '$1=="spike_in" && ($2=="43")' | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print C/(C+T)}'`
    spikein43_depth=`cat $f | awk '$1=="spike_in" && ($2=="43")' | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print (C+T)}'`
    fpratio=`cat $f |  grep 2kb_spikein | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print C/(C+T)}'`
    fpdepth=`cat $f |  grep 2kb_spikein | awk 'BEGIN{OFS="\t"}{C+=$6;T+=$5}END{print (C+T)}'`

    echo -e "$f\t${spikein29_ratio}\t${spikein29_depth}\t${spikein43_ratio}\t${spikein43_depth}\t${fpratio}\t${fpdepth}" >> all_conversion_fp.txt
done

#####
echo -e "smp\tNNUNN1_ratio\tNNUNN1_depth" >>all_conversion_fp.txt
for i in `ls H*_NNUNN1.contex.table`
do
        num_columns=$(awk '{print NF; exit}' $i)
        if [ $num_columns -eq 8 ]; then
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$7}END{print C/(C+T), (C+T)}'`
        else
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$6}END{print C/(C+T), (C+T)}'`
        fi
        echo -e "$i\t${nnunn}" >>all_conversion_fp.txt
done

echo -e "smp\tNNUNN2_ratio\tNNUNN2_depth" >>all_conversion_fp.txt
for i in `ls *_NNUNN2.contex.table`
do
        num_columns=$(awk '{print NF; exit}' $i)
        if [ $num_columns -eq 8 ]; then
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$7}END{print C/(C+T), (C+T)}'`
        else
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$6}END{print C/(C+T), (C+T)}'`
        fi
    echo -e "$i\t${nnunn}" >>all_conversion_fp.txt
done

###############
module load BEDTools/2.30.0-GCC-10.2.0
#####treat1 VS input1
zcat mergebam/treat1.mRNA_pseusite.mpile.txt.gz| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}' > mergebam/treat1.bed
bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed mergebam/treat1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> mergebam/treat1.fasta.bed
zcat mergebam/treat1.mRNA_pseusite.mpile.txt.gz | bedtools intersect -a - -b mergebam/treat1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}' > mergebam/treat1.sta.txt

bedtools intersect -a <(cat mergebam/treat1.sta.txt |sort -k1,1 -k2,2n) -b <(zcat mergebam/input1.mRNA_pseusite.mpile.txt.gz| sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$9/$6,$10,$15,$16,$17,$18,$18/($17+$18),$19,$19/$16}' > mergebam/treat1_input1_allT_motif.bed

rm mergebam/treat1.bed
rm mergebam/treat1.fasta.bed
rm mergebam/treat1.sta.txt
rm mergebam/input1.bed
#colnames chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio

####treat2 vs input2
zcat mergebam/treat2.mRNA_pseusite.mpile.txt.gz| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-2,$3+2,$4,$5,$6,$7,$8,$9}' > mergebam/treat2.bed
bedtools getfasta -fi /users/ludwig/ebu571/ebu571/resource/GRCh38.fa -bed mergebam/treat2.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> mergebam/treat2.fasta.bed
zcat mergebam/treat2.mRNA_pseusite.mpile.txt.gz | bedtools intersect -a - -b mergebam/treat2.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$9,$13}' > mergebam/treat2.sta.txt

bedtools intersect -a <(cat mergebam/treat2.sta.txt |sort -k1,1 -k2,2n) -b <(zcat mergebam/input2.mRNA_pseusite.mpile.txt.gz| sort -k1,1 -k2,2n) -sorted -wa -wb -loj | awk '($7+$8)>0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$8,$8/($7+$8),$9,$9/$6,$10,$15,$16,$17,$18,$18/($17+$18),$19,$19/$16}' > mergebam/treat2_input2_allT_motif.bed

rm mergebam/treat2.bed
rm mergebam/treat2.fasta.bed
rm mergebam/treat2.sta.txt
rm mergebam/input2.bed
















##/well/ludwig/users/ebu571/conda/skylake/envs/variant_effect/bin/variant-effect -i merge1_pseusites.tsv -H -r human -t RNA  > merge1_pseusites_annotated.tsv

#####################fing all potential snp site in input
####mutation ratio > 1% and number of nutation >=2
cat input1.mRNA_pseusite.mpile.txt | awk '($7+$8)>0'|awk '{OFS="\t"}{print $0, $8/($7+$8)}' | awk '$NF>=0.01 && $8>=2' | awk '($7+$8)>=20' | gzip > input1.background.txt.gz
cat input2.mRNA_pseusite.mpile.txt | awk '($7+$8)>0'|awk '{OFS="\t"}{print $0, $8/($7+$8)}' | awk '$NF>=0.01 && $8>=2' | awk '($7+$8)>=20' | gzip > input2.background.txt.gz

###very high confident snp mutation ratio higher than 25%
zcat input1.background.txt.gz | awk '$NF>=0.25' | wc -l ###4698
zcat input2.background.txt.gz | awk '$NF>=0.25' | wc -l ###4625

bedtools intersect -a <(zcat input1.background.txt.gz) -b <(zcat input2.background.txt.gz ) | wc -l   ###5943
bedtools intersect -a <(zcat input1.background.txt.gz | awk '$NF>=0.25') -b <(zcat input2.background.txt.gz | awk '$NF>=0.25') | wc -l ####4206


cat treat1_input1_pvalue* | grep -v ^motif | awk '{OFS="\t"}{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' | gzip > treat1_input1_tmp_site.txt.gz
cat treat2_input2_pvalue* | grep -v ^motif | awk '{OFS="\t"}{print $2,$3,$4,$5,$6,$1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' | gzip > treat2_input2_tmp_site.txt.gz


bedtools intersect -a <(zcat treat1_input1_tmp_site.txt.gz | sort -k1,1 -k2,2n) -b <(zcat treat2_input2_tmp_site.txt.gz| sort -k1,1 -k2,2n) -sorted -wa -wb | gzip > treat1_treat2_merge.txt.gz
zcat treat1_treat2_merge.txt.gz | wc -l ###1536164


zcat treat1_treat2_merge.txt.gz | awk '$22>=0.05' | awk '$17<0.01 && $16<=2' | wc -l  ###52424
zcat treat1_treat2_merge.txt.gz | awk '$46>=0.05' | awk '$41<0.01 && $40<=2' | wc -l  ###54741

zcat treat1_treat2_merge.txt.gz | awk '$22>=0.05' | awk '$17<0.01 && $16<=2' > treat1_input1_tmp_site.txt
zcat treat1_treat2_merge.txt.gz | awk '$46>=0.05' | awk '$41<0.01 && $40<=2' > treat2_input2_tmp_site.txt

bedtools intersect -a <(cat treat1_input1_tmp_site.txt |sort -k1,1 -k2,2n) -b <(cat treat2_input2_tmp_site.txt | sort -k1,1 -k2,2n) -sorted -wa -wb > merge_site.txt
cat merge_site.txt | wc -l ###15388

#########plot overlap between two replicates
#######different modification level cutoff. different percentage of overlap



echo -e "chr\tpos\tstrand\tref\talt" > merge_pseusites.tsv
cat merge_site.txt | awk '{OFS="\t"}{print $1, $3, $4, "T", "N"}'  >> merge_pseusites.tsv

/well/ludwig/users/ebu571/conda/skylake/envs/variant_effect/bin/variant-effect -i merge_pseusites.tsv -H -r human -t RNA  > merge_pseusites_annotated.tsv


