#!/bin/bash

#SBATCH --job-name test     # Name for your job
#SBATCH --output mimseq.out    # Standard out goes to this file
#SBATCH --error mimseq.err     # Standard err goes to this file
#SBATCH -p short
#SBATCH --cpus-per-task=4

module load Anaconda3/2021.05
source /apps/eb/2020b/skylake/software/Anaconda3/2021.05/etc/profile.d/conda.sh
conda activate mimseq1
mimseq --species Hsap --cluster-id 0.97 --threads 15 --max-mismatches 0.075 --control-condition control -n mim_test --out-dir mim_test_controlvstreat --max-multi 4 --remap --remap-mismatches 0.05 fastq/sample.txt


##test
##module load Anaconda3/2021.05
#conda activate 
#conda activate mimseq1
#/well/ludwig/users/ebu571/conda/skylake/envs/mimseq1/bin/mimseq --species Hsap --cluster-id 0.97 --threads 15 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir hg38_HEK239vsK562 --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt


#cat mim_test_controlvstreat/counts/Isodecoder_counts_raw.txt | cut -f1,3,4 > mim_test_controlvstreat/counts/Isodecoder_counts_raw_lowexpressed.txt

join -1 4 -2 1 <(cat mim_test_controlvstreat/counts/Isodecoder_counts_raw.txt | cut -f1,3,4,8 | sort -k 4) <(cat mim_test_controlvstreat/annotation/mim_test_isodecoderTranscripts.fa.fai | cut -f1,2 | sort -k 1) | sed 's/ /\t/g'> featureCounts/input_feature.txt

####proceed to tRNA.r script for Normalized TPM value
java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar NormalizeFasta I=tRNA_high_confident.fa O=tRNA_high_confident_tmp.fa
join -1 2 -2 2 -a1 -a2 <(cat mim_test_controlvstreat/annotation/mim_test_isodecoderTranscripts.fa | paste - - | sort -k 2) <(cat resource/tRNA_high_confident_tmp.fa | paste - - | sort -k 2) > resource/tRNA_reference_comparision.txt