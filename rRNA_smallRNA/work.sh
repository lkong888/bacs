module load samtools/1.8-gcc5.4.0
module load Bowtie2/2.4.4-GCC-10.3.0

bowtie2-build --threads 2 nnunn.fa nnunn
bowtie2-build --threads 2 spikein_rRNA.fa spikein_rRNA
bowtie2-build --threads 2 tRNA_snorna_cdhit.fa tRNA_snorna_cdhit


###check conversion and false positives
echo -e "smp\ttotal_count\tA\tG\tC\tT\tUtoC\tUtoR\tconversion\tfp_ratio\tfp_depth"> align/spike_in.txt
for i in `ls align/*.rRNA.filter.dedup.mpile.txt`
do
num=`cat $i |grep spike_in | awk '$2==29 || $2==43' | awk 'BEGIN{OFS="\t"}{A+=($5+$14);G+=($7+$16); C+=($6+$15);T+=($8+$17);gap+=($11+$20)}END{print (A+G+C+T), A, G, C, T, G/(A+G+C+T), (C+T)/(A+G+C+T)}'`
conversion=`cat $i |grep spike_in | awk '$2==29 || $2==43' | awk 'BEGIN{OFS="\t"}{A+=($5+$14);G+=($7+$16)}END{print G/(G+A)}'`
fpratio=`cat $i |  grep 2kb_spikein | awk '$3=="T"' | awk 'BEGIN{OFS="\t"}{C+=($6+$15);T+=($8+$17)}END{print C/(C+T)}'`
fpdepth=`cat $i |  grep 2kb_spikein | awk '$3=="T"' | awk 'BEGIN{OFS="\t"}{C+=($6+$15);T+=($8+$17)}END{print (C+T)}'`
echo -e "$i\t${num}\t${conversion}\t${fpratio}\t${fpdepth}">>align/spike_in.txt
done

echo -e "smp\tNNUNN1_ratio\tNNUNN1_depth" >align/nnunn1.txt
for i in `ls align/H*_NNUNN1.contex.table`
do
        num_columns=$(awk '{print NF; exit}' $i)
        if [ $num_columns -eq 8 ]; then
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$7}END{print C/(C+T), (C+T)}'`
        else
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$6}END{print C/(C+T), (C+T)}'`
        fi
        echo -e "$i\t${nnunn}" >>align/nnunn1.txt
done

#####GUUCN 
echo -e "smp\tGUUCN_ratio\tGUUCN_depth" >align/GUUCN.txt
for i in `ls align/H*_NNUNN1.contex.table`
do
        num_columns=$(awk '{print NF; exit}' $i)
        if [ $num_columns -eq 8 ]; then
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '$1=="GT"'|awk '$2=="CA" ||$2=="CT"|| $2=="CC" || $2=="CG"' |awk '{OFS="\t"}{C+=$4;T+=$7}END{print C/(C+T), (C+T)}'`
        else
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '$1=="GT"'|awk '$2=="CA" ||$2=="CT"|| $2=="CC" || $2=="CG"' |awk '{OFS="\t"}{C+=$4;T+=$6}END{print C/(C+T), (C+T)}'`
        fi
        echo -e "$i\t${nnunn}" >>align/GUUCN.txt
done



echo -e "smp\tNNUNN2_ratio\tNNUNN2_depth" >align/nnunn2.txt
for i in `ls align/*_NNUNN2.contex.table`
do
        num_columns=$(awk '{print NF; exit}' $i)
        if [ $num_columns -eq 8 ]; then
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$7}END{print C/(C+T), (C+T)}'`
        else
        nnunn=`cat $i | sed 's/ /\t/g' | grep -v ^first | awk '{OFS="\t"}{C+=$4;T+=$6}END{print C/(C+T), (C+T)}'`
        fi
    echo -e "$i\t${nnunn}" >>align/nnunn2.txt
done


##########how many psi sites on 18S, 28S and mtrRNAs? H7 vs H9; H8vs H10; H3vs H16; H23 vs H26
####18S

module load BEDTools/2.30.0-GCC-10.2.0
#####h7 VS h9
cat align/H7_2023May19_S7.rRNA.filter.dedup_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/H7_2023May19_S7.bed
bedtools getfasta -fi resource/spikein_rRNA.fa -bed align/H7_2023May19_S7.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/H7_2023May19_S7.fasta.bed
cat align/H7_2023May19_S7.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/H7_2023May19_S7.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/H7_2023May19_S7.sta.txt

cat align/H9_2023May19_S9.rRNA.filter.dedup_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/H9_2023May19_S9.bed
bedtools intersect -a align/H7_2023May19_S7.sta.txt -b align/H9_2023May19_S9.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/H7_H9_allT_motif.bed

#colnames chr, start, end, ref, depth, T, C, mod_level, gap, gap_ratio, motif, ctrl_depth, ctrl_T, ctrl_C, ctrl_ratio, ctrl_gap, ctrl_gap_ratio

####H8 vs H10
cat align/H8_2023May19_S8.rRNA.filter.dedup_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/H8_2023May19_S8.bed
bedtools getfasta -fi resource/spikein_rRNA.fa -bed align/H8_2023May19_S8.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/H8_2023May19_S8.fasta.bed
cat align/H8_2023May19_S8.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/H8_2023May19_S8.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/H8_2023May19_S8.sta.txt

cat align/H10_2023May19_S10.rRNA.filter.dedup_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/H10_2023May19_S10.bed
bedtools intersect -a align/H8_2023May19_S8.sta.txt -b align/H10_2023May19_S10.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/H8_H10_allT_motif.bed


####H3 vs H116
cat align/H3_2023Apr16_S3.rRNA.filter.dedup_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/H3_2023Apr16_S3.bed
bedtools getfasta -fi resource/spikein_rRNA.fa -bed align/H3_2023Apr16_S3.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/H3_2023Apr16_S3.fasta.bed
cat align/H3_2023Apr16_S3.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/H3_2023Apr16_S3.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/H3_2023Apr16_S3.sta.txt

cat align/H16_2023Apr16_S15.rRNA.filter.dedup_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/H16_2023Apr16_S15.bed
bedtools intersect -a align/H3_2023Apr16_S3.sta.txt -b align/H16_2023Apr16_S15.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/H3_H16_allT_motif.bed


####H23 vs H26
cat align/H23_2023Apr16_S19.rRNA.filter.dedup_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/H23_2023Apr16_S19.bed
bedtools getfasta -fi resource/spikein_rRNA.fa -bed align/H23_2023Apr16_S19.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/H23_2023Apr16_S19.fasta.bed
cat align/H23_2023Apr16_S19.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/H23_2023Apr16_S19.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/H23_2023Apr16_S19.sta.txt

cat align/H26_2023Apr16_S22.rRNA.filter.dedup_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/H26_2023Apr16_S22.bed
bedtools intersect -a align/H23_2023Apr16_S19.sta.txt -b align/H26_2023Apr16_S22.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/H23_H26_allT_motif.bed


###for snoRNA and tRNA sites

#treat1 vs input1
cat align/treat1.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/treat1.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed align/treat1.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/treat1.fasta.bed
cat align/treat1.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/treat1.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/treat1.sta.txt

cat align/input1.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/input1.bed
bedtools intersect -a align/treat1.sta.txt -b align/input1.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/treat1_input1_motif.bed

#treat2 vs input2
cat align/treat2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/treat2.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed align/treat2.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/treat2.fasta.bed
cat align/treat2.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/treat2.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/treat2.sta.txt

cat align/input2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/input2.bed
bedtools intersect -a align/treat2.sta.txt -b align/input2.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/treat2_input2_motif.bed

#treat3 vs input2
cat align/treat3.snoRNA.tRNA.filter.sort_pseusite.mpile.txt| awk '$2>3' | awk '{OFS="\t"}{print $1, $2-3,$2+2, $3,$4,$5,$6,$7}' |sed 's/-/_/g' > align/treat3.bed
bedtools getfasta -fi resource/tRNA_snorna.fa -bed align/treat3.bed | paste - - | sed 's/>//g' | sed 's/:/\t/g'| sed 's/-/\t/g' | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4}'> align/treat3.fasta.bed
cat align/treat3.bed | awk '{OFS="\t"}{print $1, $2+2, $3-2, $4,$5,$6,$7,$8}' | bedtools intersect -a - -b align/treat3.fasta.bed -wa -wb -loj | \
        awk '{OFS="\t"}{print $1, $2, $3, $4,$5,$6,$7,$8,$12}' > align/treat3.sta.txt

cat align/input2.snoRNA.tRNA.filter.sort_pseusite.mpile.txt | sed 's/-/_/g' | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $5, $6,$7}' > align/input2.bed
bedtools intersect -a align/treat3.sta.txt -b align/input2.bed -wa -wb -loj | awk '($6+$7)>0' | awk '($15+$16) >0'| awk '{OFS="\t"}{print $1, $2,$3,$4,$5,$6,$7,$7/($6+$7),$8,$8/$5,$9,$14,$15,$16, $16/($15+$16),$17,$17/$14}' > align/treat3_input2_motif.bed



##########exclude some anticodons
join -1 1 -2 1 -v1 <( cat snoRNA_tRNA_called_sites_fdr_filtered.txt | sort -k1 ) <( cat input2_TPM_lower_quantile.txt | sort -k1 ) | sed 's/ /\t/g' > snoRNA_tRNA_called_sites_fdr_filtered_tpmexclude.txt ###635



#################snoRNA annotation###########################
cat /users/ludwig/ebu571/ebu571/resource/snoRna_target.txt /users/ludwig/ebu571/ebu571/resource/snoRna_target1.txt | sed 's/,/\t/g' | awk '$4~"NR"' | cut -f2,4 > resource/snoRNA_index.txt
join -1 1 -2 1 <(cat resource/snoRNA_index.txt | sort -k1 ) \
        <(cat /users/ludwig/ebu571/ebu571/resource/snodb.annotation.txt /users/ludwig/ebu571/ebu571/resource/snodb.annotation1.txt | sed 's/,/\t/g' | sort -k1 ) > resource/snoRNA_annptation.txt

module load BEDTools/2.30.0-GCC-10.2.0
cat site/snosite_namechanged.txt | grep -v ^chr |awk '{OFS="\t"}{print $NF, $2,$3, $4,$5,$6,$7,$8, $9,$1}'
bedtools map -a <(cat resource/snoRNA_annptation.txt | sed 's/ /\t/g' | awk '{OFS="\t"}{print $2, $3, $4+1, $5, $6, $1}'| sort -k1,1 -k2,2n) -b <(cat site/snosite_namechanged.txt | grep -v ^chr |awk '{OFS="\t"}{print $NF, $2,$3, $4,$5,$6,$7,$8, $9,$1}'| sort -k1,1 -k2,2n) -c 9 -o mean 

bedtools intersect -a <(cat resource/snoRNA_annptation.txt | sed 's/ /\t/g' | awk '{OFS="\t"}{print $2, $3, $4+1, $5, $6, $1}'| sort -k1,1 -k2,2n) -b <(cat site/snosite_namechanged.txt | grep -v ^chr |awk '{OFS="\t"}{print $NF, $2,$3, $4,$5,$6,$7,$8, $9,$1}'| sort -k1,1 -k2,2n) -wa -wb > site/sno_annotated_site.txt
 


######for snoRNA, divide each reference into 100 bins
bedtools makewindows -g resource/snoRNA_reference.size.txt -n 40 -i winnum > resource/snoRNA_40bins.txt
bedtools map -a <(cat resource/snoRNA_40bins.txt |sort -k1,1 -k2,2n) -b <(cat site/snoRNA_tRNA_called_sites_fdrpassed.txt | grep -v tRNA | grep -v ^chr| sort -k1,1 -k2,2n) -o count > site/snoRNA_density_40bins.txt       #chr, start,end,binnum,count