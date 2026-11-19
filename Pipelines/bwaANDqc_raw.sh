## define variables
fasta=/data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
readdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/fastq
mapdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/bam
samples=$(ls  $readdir/*.fastq)
MINMAPQ=30

## Mapping with BWA 
{ # build bwa index
bwa index /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
cd /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/

# mapping
for fq in $samples;
do
    mapdir_sample=$(basename $fq)
    mapdir_sample=${mapdir}/${mapdir_sample/.fastq/}
   
  OUTPREF=$mapdir_sample
  
  OUTSAM=$(mktemp -p $(dirname ${OUTPREF}))
  SAI=$(mktemp -p $(dirname ${OUTPREF}))
    # Single End mode
    echo "single-end mode"
    echo $fq
    /cellfile/datapublic/tpadvits/SOFT/bwa-0.7.17/bwa aln -t 48 $fasta $fq > $SAI
    /cellfile/datapublic/tpadvits/SOFT/bwa-0.7.17/bwa samse $fasta $SAI $fq > $OUTSAM
   
    #echo $CMD
    eval $CMD
  
  ## Data conversion
  samtools view --threads 48 -S -b $OUTSAM > ${OUTPREF}_raw.bam
  samtools sort --threads 48 -o ${OUTPREF}_.bam ${OUTPREF}_raw.bam 
  samtools index ${OUTPREF}_.bam
  samtools flagstat ${OUTPREF}_.bam > ${OUTPREF}.stats
  ## Data filtering
  # remove poor MAPQ, then remove duplicates
  samtools view --threads 48 -b -q $MINMAPQ ${OUTPREF}_.bam > ${OUTPREF}_MAPQ30.bam
  samtools index ${OUTPREF}_MAPQ30.bam 
  samtools rmdup -s ${OUTPREF}_MAPQ30.bam ${OUTPREF}_MAPQ30unique.bam
  samtools index ${OUTPREF}_MAPQ30unique.bam
  samtools flagstat ${OUTPREF}_MAPQ30unique.bam > ${OUTPREF}_MAPQ30unique.stats 

  rm $SAI $OUTSAM ${OUTPREF}_raw.bam
  #echo ${OUTPREF}_filtered.bam
done
}


# # plot MAPQ

# for F in $(ls $mapdir/*_.bam); 
# do
# 
#   BAM=$F
#   OUTPREF=$mapdir/images/$(basename $F)
#  
#   OUTDIR=$(dirname $OUTPREF)
#   if [[ ! -d $OUTDIR ]];
#   then
#     echo -e "$OUTDIR missing. \nExit..."
#     exit
#   fi
#   TMP=$(mktemp -p $(dirname $OUTPREF) )
# 
#   echo "Creating MAPQ histogram"
#   samtools view $BAM | cut -f 5 > $TMP
# #   
#   Rscript -e "args=commandArgs(T); out=paste(args[2],'_distr.png', sep =''); png(out, width = 1024, height=768, pointsize=18); x=as.matrix(read.table(args[1])) ;hist(x, main = paste(args[3], sep =''), xlab = 'MAPQ', breaks=seq(0,max(x), 1)); dev.off();" $TMP $OUTPREF $(basename $BAM)
#   
#   rm $TMP
# done
#}


# ### Deeptool QC
#{ 
# bam=$(ls  $mapdir/*_.bam)
# bam30=$(ls  $mapdir/*MAPQ30unique.bam)
#    
# # /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/multiBamSummary bins --bamfiles $bam --labels WT1martin_IPrep1 WT1martin_IPrep2 WT1martin_IPrep3 WT1martin_input_rep1 WT1martin_input_rep2 WT1martin_input_rep3 -out $mapdir/deeptools/mbamsum.npz 
# /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/multiBamSummary bins --bamfiles $bam30 --labels WT1martin_IPrep1 WT1martin_IPrep2 WT1martin_IPrep3 WT1martin_input_rep1 WT1martin_input_rep2 WT1martin_input_rep3 -out $mapdir/deeptools/MQ30uni_mbamsum.npz 
# 
# # /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotCorrelation \
# #     -in $mapdir/deeptools/mbamsum.npz \
# #     --corMethod spearman --skipZeros \
# #     --plotTitle "Spearman Correlation of Read Counts" \
# #     --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
# #     -o $mapdir/deeptools/heatmap_SpearmanCorr_readCounts.png \
# #     --outFileCorMatrix $mapdir/deeptools/SpearmanCorr_readCounts.tab  
#     
# /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotCorrelation \
#     -in $mapdir/deeptools/MQ30uni_mbamsum.npz \
#     --corMethod spearman --skipZeros \
#     --plotTitle "Spearman Correlation of Read Counts, MAPQ30uni" \
#     --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
#     -o $mapdir/deeptools/heatmap_SpearmanCorr_MQ30uni_readCounts.png \
#     --outFileCorMatrix $mapdir/deeptools/SpearmanCorr_MQ30uni_readCounts.tab    
#     
# # /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotPCA -in $mapdir/deeptools/mbamsum.npz \
# # -o $mapdir/deeptools/ChipseqPCA_readCounts.png \
# # -T "PCA of read counts"
# 
# /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotPCA -in $mapdir/deeptools/MQ30uni_mbamsum.npz \
# -o $mapdir/deeptools/ChipseqPCA_MQ30uni_readCounts.png \
# -T "PCA of read counts, MAPQ30uni"
# 
# /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotFingerprint -b $bam \
# --labels IPrep1 IPrep2 IPrep3 input_rep1 input_rep2 input_rep3 \
# --skipZeros \
# -T "Fingerprints of Chip-seq samples"  \
# --plotFile $mapdir/deeptools/Chipseq_fingerprints.png
# 
# /cellfile/datapublic/tpadvits/SOFT/deepTools/bin/plotFingerprint -b $bam30 \
# --labels IPrep1 IPrep2 IPrep3 input_rep1 input_rep2 input_rep3 \
# --skipZeros \
# -T "Fingerprints of Chip-seq samples, MAPQ30"  \
# --plotFile $mapdir/deeptools/Chipseq_fingerprintsMAPQ30.png
#}

# ### cross correlation analysis with spp #########
# { 
# mapdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30
# outdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30/spp_crosscorr
# mkdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/bam/spp_crosscorr
# 
# bam_filt=$(ls  $mapdir/*filtSE.bam)
# for bm in $bam_filt;
#      do    
# Rscript /data/user/tpadvits/BIOSOFT/ccQualityControl.v.1.1/phantompeakqualtools/run_spp_nodups.R -c=$bm -savp -odir=$outdir -p=8 -out=spp_crosscorr.txt
#      done
# }      

 ### Call peaks with MACS2 
#{
#peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/peaks/MACS2

# NAME="WT1martinSE_rep1_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703423_MAPQ30unique.bam -c $mapdir/SRR1703426_MAPQ30unique.bam -n $NAME -B -g mm -q 0.05 -s 40 --call-summits --outdir $peakdir
#     
# NAME="WT1martin_SE_rep2_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703424_MAPQ30unique.bam -c $mapdir/SRR1703427_MAPQ30unique.bam -n $NAME -B -g mm -q 0.05 -s 40 --call-summits --outdir $peakdir
# 
# NAME="WT1martin_SE_rep3_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703425_MAPQ30unique.bam -c $mapdir/SRR1703428_MAPQ30unique.bam -n $NAME -B -g mm -q 0.05 -s 40 --call-summits --outdir $peakdir
# 
#     
# ## call peaks for pooled IP/input
# samtools merge $mapdir/WT1martin_merged_Input.bam $mapdir/SRR1703426_MAPQ30unique.bam $mapdir/SRR1703427_MAPQ30unique.bam $mapdir/SRR1703428_MAPQ30unique.bam 
# samtools merge $mapdir/WT1martin_merged_IP.bam  $mapdir/SRR1703423_MAPQ30unique.bam  $mapdir/SRR1703424_MAPQ30unique.bam $mapdir/SRR1703425_MAPQ30unique.bam 
# 
# peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/peaks/MACS2/merged
# mkdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/peaks/MACS2/merged
# 
# NAME="WT1martin_SE_pooled_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/WT1martin_merged_IP.bam -c $mapdir/WT1martin_merged_Input.bam -n $NAME -B -g mm -q 0.05 -s 40 --keep-dup all --call-summits --outdir $peakdir
# 

# 
# ####### do IDR analysis #######
# ### Call low confedence peaks with MACS2 
# peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/peaks/MACS2/p0.01
# mkdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/peaks/MACS2/p0.01
# 
# NAME="WT1martinSE_rep1p0.01_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703423_MAPQ30unique.bam -c $mapdir/SRR1703426_MAPQ30unique.bam -n $NAME -B -g mm -p 0.01 -s 40 --call-summits --outdir $peakdir
#     
# NAME="WT1martin_SE_rep2p0.01_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703424_MAPQ30unique.bam -c $mapdir/SRR1703427_MAPQ30unique.bam -n $NAME -B -g mm -p 0.01 -s 40 --call-summits --outdir $peakdir
# 
# NAME="WT1martin_SE_rep3p0.01_14.05.18"
#     macs2 callpeak -f BAM -t $mapdir/SRR1703425_MAPQ30unique.bam -c $mapdir/SRR1703428_MAPQ30unique.bam -n $NAME -B -g mm -p 0.01 -s 40 --call-summits --outdir $peakdir


# # do IDR (locally)
# 
# PATH=$PATH:~/BIOSOFT/data/user/tpadvits/BIOSOFT/idr-2.0.2/bin
# module unload python
# module unload python3
# module unload python-p
# module load python3
# /cellfile/datapublic/tpadvits/SOFT/idr-2.0.2/bin/idr --samples $peakdir/WT1martin_SE_rep1p0.01_14.05.18_peaks.narrowPeak $peakdir/WT1martin_SE_rep2p0.01_14.05.18_peaks.narrowPeak  --input-file-type narrowPeak --plot -o $peakdir/MACSp0.01_rep12_idrValues
# /cellfile/datapublic/tpadvits/SOFT/idr-2.0.2/bin/idr --samples $peakdir/WT1martin_SE_rep2p0.01_14.05.18_peaks.narrowPeak $peakdir/WT1martin_SE_rep3p0.01_14.05.18_peaks.narrowPeak  --input-file-type narrowPeak --plot -o $peakdir/MACSp0.01_rep23_idrValues
# /cellfile/datapublic/tpadvits/SOFT/idr-2.0.2/bin/idr --samples $peakdir/WT1martin_SE_rep1p0.01_14.05.18_peaks.narrowPeak $peakdir/WT1martin_SE_rep3p0.01_14.05.18_peaks.narrowPeak  --input-file-type narrowPeak --plot -o $peakdir/MACSp0.01_rep13_idrValues
#}

# ### subsample and call peaks on subsampled pseudo replicates
#{
# module unload python3
# module load python
# 
# bam30=$(ls  $mapdir/*MAPQ30unique.bam)
# IP=$(echo $bam30 | cut -d" " -f -3)
# for F in $IP;
# do
# out=$(basename $F)
# out=${out/_MAPQ30unique.bam/}
# len=$(samtools view -c $F)
# len=$(echo $((len / 2)))
# macs2 randsample -t $F -n $len -o $mapdir/${out}_pr1.bam
# macs2 randsample -t $F -n $len -o $mapdir/${out}_pr2.bam
# done
# 
# merge pseudoreplicates in pooled pseudoreplicates
# cat $mapdir/SRR1703423_pr1.bam $mapdir/SRR1703424_pr1.bam $mapdir/SRR1703425_pr1.bam > $mapdir/pooledIP_pseudo_rep1.bed
# cat $mapdir/SRR1703423_pr2.bam $mapdir/SRR1703424_pr2.bam $mapdir/SRR1703425_pr2.bam > $mapdir/pooledIP_pseudo_rep2.bed
# 
# convert input from bam to bed
# bedtools bamtobed -i $mapdir/WT1martin_merged_Input.bam >  $mapdir/WT1martin_merged_Input.bed
# 
# call peaks on pooled pseudo replicates
# peakdir=$mapdir/pseudo_rep
# mkdir $peakdir
# module unload python3
# module load python
# NAME="WT1martinSE_ppr1p0.01_14.05.18"
# macs2 callpeak -f BED -t $mapdir/pooledIP_pseudo_rep1.bed -c $mapdir/WT1martin_merged_Input.bed -n $NAME -B -g mm -p 0.01 -s 40 --keep-dup all --call-summits --outdir $peakdir
# NAME="WT1martinSE_ppr2p0.01_14.05.18"
# macs2 callpeak -f BED -t $mapdir/pooledIP_pseudo_rep2.bed -c $mapdir/WT1martin_merged_Input.bed -n $NAME -B -g mm -p 0.01 -s 40 --keep-dup all --call-summits --outdir $peakdir
# 
# do IDR on peaks called from pooled pseudo replicates
# module unload python
# module unload python3
# module unload python-p
# module load python3
# /cellfile/datapublic/tpadvits/SOFT/idr-2.0.2/bin/idr --samples $peakdir/WT1martinSE_ppr1p0.01_14.05.18_peaks.narrowPeak $peakdir/WT1martinSE_ppr2p0.01_14.05.18_peaks.narrowPeak  --input-file-type narrowPeak --plot -o $peakdir/MACSp0.01_ppr_idrValues
#}


# ###### HOMER motif analysis ##########
# { 
# PATH=$PATH:/data/public/tpadvits/SOFT/homer.v4.9.1/bin
# 
# mot_dir=/data/public/tpadvits/global_data/motif_analysis
# 
# # /cellfile/datapublic/tpadvits/SOFT/homer.v4.9.1/bin/findMotifsGenome.pl $mapdir/pseudo_rep/MACSp0.01_ppr_idr0.05_UCSC.csv mm10 $mapdir/pseudo_rep/HOMER/WT1martin_HOCOMOCO_HOMER2 -size 200 -len 8,10,12,14 -p 32 -mcheck $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif -mknown $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif 
# 
# 
# dir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Martin2014_collaborator/Data/qc/logs
# 
# # findMotifsGenome.pl $dir/WT1martin_tead1_both_ench.csv mm10 $dir/WT1martin_tead1_ench_HOCOMOCO_HOMER2 -size 200 -len 8,10,12,14 -p 32 -mcheck $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif -mknown $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif 
# 
# findMotifsGenome.pl $dir/WT1martin_tead1_both_prom.csv mm10 $dir/WT1martin_tead1_prom_HOCOMOCO_HOMER2 -size 200 -len 8,10,12,14 -p 32 -mcheck $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif -mknown $mot_dir/HOCOMOCOv11mouse_HOMER_all.motif 
# }


#========= browser plot ===============
# prepare data for plotting  https://github.com/taoliu/MACS/wiki/Build-Signal-Track 
mapdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30
# samtools merge -b $mapdir/SRR1703423_pr1.bam $mapdir/SRR1703424_pr1.bam $mapdir/SRR1703425_pr1.bam > $mapdir/pooledIP.bam
# cat $mapdir/SRR1703423_pr2.bam $mapdir/SRR1703424_pr2.bam $mapdir/SRR1703425_pr2.bam > $mapdir/pooledIP_pseudo_rep2.bed


## difference (Fold Enrichment) between 2 bedgraphs using macs2 and separate chip-seq replicates 
## Wt1wt
# run Macs to to generate bedgraph files, no need to pool bam files! 
macs2 callpeak -t $mapdir/SN7640223_16878_ChipWt1wtF0513Rep3_MAPQ30filtSE.bam $mapdir/SN7640223_16880_ChipWt1wtI2013Rep1_MAPQ30filtSE.bam -c $mapdir/SN7640223_16879_InputWt1wtRep3_MAPQ30filtSE.bam $mapdir/SN7640223_16881_InputWt1wtRep1_MAPQ30filtSE.bam  -B --nomodel --extsize 147 --SPMR -g mm --keep-dup 2 -n $mapdir/Wt1wt_normRep
# bedgraph coverage normalized to input
macs2 bdgcmp -t $mapdir/Wt1wt_normRep_treat_pileup.bdg -c $mapdir/Wt1wt_normRep_control_lambda.bdg -o $mapdir/Wt1wt_IPoverInput_repFE.bdg -m FE
# bedgraph to bigwig using UCSC utility
/data/user/tpadvits/BIOSOFT/UCSCutils/bedGraphToBigWig $mapdir/Wt1wt_IPoverInput_repFE.bdg /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/chrNameLength.txt $mapdir/Wt1wt_normSignal_rep.bw

## Wt1het
macs2 callpeak -t $mapdir/SN7640223_16876_ChipWt1hetF0513Rep3_MAPQ30filtSE.bam $mapdir/SN7640285_29526_ChipWt1hetI2013Rep1_MAPQ30filtSE.bam -c   $mapdir/WT1het_input_merged.bam $mapdir/SN7640223_16877_InputWt1hetRep3_MAPQ30filtSE.bam $mapdir/SN7640285_29527_InputWt1hetRep1_MAPQ30filtSE.bam  -B --nomodel --extsize 147 --SPMR -g mm --keep-dup 2 -n  $mapdir/Wt1het_normRep
macs2 bdgcmp -t $mapdir/Wt1het_normRep_treat_pileup.bdg -c $mapdir/Wt1het_normRep_control_lambda.bdg -o $mapdir/Wt1het_IPoverInput_repFE.bdg -m FE
/data/user/tpadvits/BIOSOFT/UCSCutils/bedGraphToBigWig $mapdir/Wt1het_IPoverInput_repFE.bdg /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/chrNameLength.txt $mapdir/Wt1het_normSignal_rep.bw

