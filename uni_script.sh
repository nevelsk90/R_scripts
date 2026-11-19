PATH=$PATH:~/BIOSOFT/sratoolkit.2.8.0-ubuntu64/bin
PATH=$PATH:~/BIOSOFT/kallisto_linux-v0.43.0
PATH=$PATH:~/BIOSOFT/FastQC
PATH=$PATH:~/BIOSOFT/samstat-1.5.1/src
PATH=$PATH:~/BIOSOFT/Trimmomatic-0.36
PATH=$PATH:~/BIOSOFT/flexbar-2.7.0
PATH=$PATH:~/BIOSOFT/bowtie2-2.3.1
PATH=$PATH:~/BIOSOFT/STAR/bin/Linux_x86_64
PATH=$PATH:~/BIOSOFT/STAR/bin/Linux_x86_64
PATH=$PATH:~/BIOSOFT/bin
PATH=$PATH:~/BIOSOFT/PrimerSeq-master
PATH=$PATH:~/BIOSOFT/UCSCutils
PATH=$PATH:~/BIOSOFT/bedops_linux_x86_64-v2.4.26/bin


# to fetch DNase-seq from SRA
fastq-dump  --outdir ~/ws/SRA_downloads/files -I -W  SRR231120
# to fetch RNA-seq from SRA
fastq-dump --outdir ~/ws/SRA_downloads/files --split-3 --skip-technical -I -W  SRR4422414

# to download RNA_seq from http rep.
wget -r --no-parent --http-user=benzing --http-password="$SEQ_PASSWORD" http://bifacility.uni-koeln.de/labs/benzing/downloads/martin_may_2016/
wget -r --no-parent --http-user=beyer --http-password="$SEQ_PASSWORD" http://bifacility.uni-koeln.de/labs/beyer/downloads/Tim/
wget -r --no-parent http://127.0.0.1:42427/

# to do quality check
	# 'find' command to find all bam files and to apply samstat to them
find /data/user/tpadvits/ws/Podocyte_data/RNA_seq/martin_may_2016 -type f -name "*.bam"  -exec ./samstat {} +

# convert bam to bed (for Crunch)
	bamToBed -i 569_16878_ChipWt1wtF0513Rep3_mm9.bam > 569_16878_ChipWt1wtF0513Rep3_mm9.bed

### use of modified Michael Chip-seq pipeline
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/runPeakCalling_pipeline_human.slurm 29 Padvitski02.01_BATF_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000NSW.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000NSX.fastq.gz /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000OCS_ctl.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000OCU_ctl.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000OCW_ctl.fastq.gz /data/public/tpadvits/global_data/hg19/hg19 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/processedBATF


# to run fastQC from any location
sudo ln -s /data/user/tpadvits/BIOSOFT/FastQC /usr/local/bin/fastqc
/data/user/tpadvits/BIOSOFT/FastQC/fastqc 

# run trimmomatic
java -jar trimmomatic-0.36.jar PE -threads 30 -phred64 input_forward.fq.gz input_reverse.fq.gz -baseout  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## vew head of gzip file
  fq="ENCFF000OHO.fastq.gz"
  mv "$fq" "$fq~" &&
  gzip -cd "$fq~" | sed 's/\./N/g' | gzip > "$fq"
  
### merge BAM replicates
samtools merge /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RelB_merged.bam /data/user/tpadvits/ws/Janis_PJdata/Chip-seq/harvard2cologne/RelB_R1.bam /data/user/tpadvits/ws/Janis_PJdata/Chip-seq/harvard2cologne/RelB_R2.bam 
samtools merge /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/P52AB/p52Ab_merged.bam /data/user/tpadvits/ws/Janis_PJdata/Chip-seq/harvard2cologne/p52Ab_R1.bam /data/user/tpadvits/ws/Janis_PJdata/Chip-seq/harvard2cologne/p52Ab_R2.bam

#### run MACS2 separately
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF/mapping/ctrl_filtered.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELA/RelA_R1.bam Padvitski19.01_RELA_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELA
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF/mapping/ctrl_filtered.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELB/RelB_merged.bam Padvitski19.01_RELB_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELB
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF/mapping/ctrl_filtered.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/P52AB/p52Ab_merged.bam Padvitski19.01_P52AB_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/P52AB


fastq-dump --outdir ~/ws/Podocyte_data/RNA_seq/WT1_RNAseq --split-3 --skip-technical -I -W  SRR1703430
fastq-dump --outdir ~/ws/Podocyte_data/RNA_seq/WT1_RNAseq --split-3 --skip-technical -I -W  SRR1703429
fastq-dump --outdir ~/ws/Podocyte_data/RNA_seq/WT1_RNAseq --split-3 --skip-technical -I -W  SRR1703431

PATH=$PATH:~/BIOSOFT/kallisto_linux-v0.43.0
kallisto index --make-unique -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.introns_only.idx /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/introns/introns_small.fa

data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/WT1/glomeruli_Ribo_minus/trimmed
index=/data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.introns_only.idx
cd $data_dir
samples=$(find * -type f -name "*_1P")
outdir=/data/user/tpadvits/PROJECTS/MARTIN_PJ/RNAseq_analysis/WT1/kallisto_introns_reverse
for fq in $samples;
do
  fq2=${fq/_1P/_2P}
  dirname=${fq/filtered_1P/counted_introns_r}
  echo $dirname
  kallisto quant -i $index --rf-stranded -o $outdir/$dirname -b 100 --bias -t 7 $fq $fq2 
done


##### chip-seq analysis with IgG controls
#######################################################################  
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/runPeakCalling_pipeline_human.slurm 29 Padvitski31.01_BATF_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000NSW.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/BATF/ENCFF000NSX.fastq.gz /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/EP300/ENCFF000VWD_ctl.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/EP300/ENCFF000VWK_ctl.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RAD21/ENCFF000VWF_ctl.fastq.gz,/data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RAD21/ENCFF000VWH_ctl.fastq.gz /data/public/tpadvits/global_data/hg19/hg19 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF_IgG

sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELA/RelA_R1.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF_IgG/mapping/ctrl_filtered.bam  Padvitski31.01_RELA_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedRELA_IgG



##### chip-seq analysis with Snyder inputDNA controls
#######################################################################
sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF/mapping/sample_filtered.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedCEBPZ/mapping/ctrl_filtered.bam  Padvitski01.02_BATF_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedBATF_inputSny

sbatch /cellfile/datapublic/tpadvits/PROJECTS/Side/Janis_TF/chip-seq-peak-calling_public/bin/stepwise/runMACS2.slurm /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/RELA/RelA_R1.bam /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedCEBPZ/mapping/ctrl_filtered.bam  Padvitski01.02_RELA_wt_GM12878 /data/public/tpadvits/PROJECTS/Janis_TF/Chip_seq/processedRELA_inputSny




###########
$ curl -s  "http://rest.kegg.jp/list/organism" | grep "Mus musculus"
T01002	mmu	Mus musculus (mouse)	Eukaryotes;Animals;Vertebrates;Mammals
# get all KGML files from all pathway of M.musculus
$ curl "http://rest.kegg.jp/list/pathway/T01002" | cut -f 1 | while read A; do  curl -o "${A}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done


#### PEAR merging paired ends
####################################################################### 
cd /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/bifacility.uni-koeln.de
pear --threads 11 -f SN7640223_16876_ChipWt1hetF0513Rep3_1_sequence.fq.gz -r SN7640223_16876_ChipWt1hetF0513Rep3_2_sequence.fq.gz -o ChipWt1hetF0513Rep3_merged 
pear --threads 11 -f SN7640223_16877_InputWt1hetRep3_1_sequence.fq.gz -r SN7640223_16877_InputWt1hetRep3_2_sequence.fq.gz -o InputWt1hetRep3_merged 
pear --threads 11 -f SN7640223_16878_ChipWt1wtF0513Rep3_1_sequence.fq.gz -r SN7640223_16878_ChipWt1wtF0513Rep3_2_sequence.fq.gz -o ChipWt1wtF0513Rep3_merged
pear --threads 11 -f SN7640223_16879_InputWt1wtRep3_1_sequence.fq.gz -r SN7640223_16879_InputWt1wtRep3_2_sequence.fq.gz -o InputWt1wtRep3_merged
pear --threads 11 -f SN7640223_16880_ChipWt1wtI2013Rep1_1_sequence.fq.gz -r SN7640223_16880_ChipWt1wtI2013Rep1_2_sequence.fq.gz -o ChipWt1wtI2013Rep1_merged
pear --threads 11 -f SN7640223_16881_InputWt1wtRep1_1_sequence.fq.gz -r SN7640223_16881_InputWt1wtRep1_2_sequence.fq.gz -o InputWt1wtRep1_merged



##### c elegans vhl1 analysis
PATH=$PATH:~/BIOSOFT/kallisto_linux-v0.43.0
kallisto index -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Celegans.WBcel235.cdna_ncRNA.all.idx /data/public/tpadvits/global_data/Genome_annot/c_elegans/Caenorhabditis_elegans.WBcel235.cdna_ncRNA.all.fa


data_dir=/cellfile/datapublic/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_polyA/raw_data
cd $data_dir
samples=$(find * -type f -name "*_1_sequence*")
for fq in $samples;
do
  fq2=${fq/_1_sequence/_2_sequence}
  dirname=${fq/_1_sequence.fq.gz/_counted}
  echo $dirname
  kallisto quant -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Celegans.WBcel235.cdna_ncRNA.all.idx -o /data/public/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_ribominus/Kallisto_quant/$dirname -b 100 --bias -t 8 $fq $fq2
done

# single-end library
data_dir=/cellfile/datapublic/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_polyA/raw_data
cd $data_dir
samples=$(find * -type f)
for fq in $samples;
do
  dirname=${fq/.fastq.gz/_counted}
  echo $dirname
  kallisto quant -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Celegans.WBcel235.cdna_ncRNA.all.idx -o /data/public/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_polyA/Kallisto_quant/$dirname --single -l 96 -s 1 -b 100 --bias -t 8 $fq
done

kallisto quant -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Celegans.WBcel235.cdna_ncRNA.all.idx -o /cellfile/datapublic/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_polyA --single -l 96 -s 1 --rf-stranded -b 100 --bias -t 8 /cellfile/datapublic/tpadvits/PROJECTS/MSc/timeseries_RNAseq/timeseries_polyA/raw_data/polyA_N2day1_rep1.fastq.gz

###########################
## build bowtie indexes for mouse
bowtie2-build /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.chromosome.1.fa /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/bowtie_index/Mus_musculus.GRCm38.dna.primary


##########################################
### run pipeline_1 for mapping and QC ####
sbatch /data/public/tpadvits/PROJECTS/PodocytePJ/Chip_seq/Pipeline/pipeline_1.slurm 29  /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/bifacility.uni-koeln.de /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/bowtie_index/Mus_musculus.GRCm38.dna.primary /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/mapping_QC_results

    adapters="/data/public/tpadvits/global_data/adapters/TruSeq3-SE.fa:2:30:12"
  
java -jar /cellfile/datapublic/tpadvits/SOFT/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -trimlog /data/public/tpadvits/PROJECTS/PodocytePJ/Chip_seq/WT1/mapping_QC_result/fastq_filtered2/SN7640285_29524_ChIPWt1wtF0513Rep2_trimmed.log -$score /cellfile/DATASETS/Mouse/Kidney/Chip_seq/Bulk/In_house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/fastq/Wt1/SN7640285_29524_ChIPWt1wtF0513Rep2_1_sequence.fq.gz /data/public/tpadvits/PROJECTS/PodocytePJ/Chip_seq/WT1/mapping_QC_result/fastq_filtered2/SN7640285_29524_ChIPWt1wtF0513Rep2_filtered ILLUMINACLIP:$adapters $options

 
peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/peaks/Wt1/peaks_SEbwa_MAPQ30
bamdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30
cd $peakdir

peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/peaks/Wt1/peaks_SEbwa_MAPQ30

NAME="WT1_wt_SE_rep2p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640285_29524_ChIPWt1wtF0513Rep2_MAPQ30filtSE.bam -c $bamdir/SN7640285_29525_InputWt1wtRep2_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
NAME="WT1_het_SE_rep2p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640285_29522_ChIPWt1hetF0513Rep2_MAPQ30filtSE.bam -c $bamdir/SN7640285_29523_InputWt1hetRep2_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
    
NAME="WT1_het_SE_rep3p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640223_16876_ChipWt1hetF0513Rep3_MAPQ30filtSE.bam -c $bamdir/SN7640223_16877_InputWt1hetRep3_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
    
NAME="WT1_wt_SE_rep3p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640223_16878_ChipWt1wtF0513Rep3_MAPQ30filtSE.bam -c $bamdir/SN7640223_16879_InputWt1wtRep3_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
 
NAME="WT1_wt_SE_rep1p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640223_16880_ChipWt1wtI2013Rep1_MAPQ30filtSE.bam -c $bamdir/SN7640223_16881_InputWt1wtRep1_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
    
NAME="WT1_het_SE_rep1p0.001_04.01.18"
    macs2 callpeak -f BAM -t $bamdir/SN7640285_29526_ChipWt1hetI2013Rep1_MAPQ30filtSE.bam -c $bamdir/SN7640285_29527_InputWt1hetRep1_MAPQ30filtSE.bam -n $NAME -B -g mm -p 0.001 --call-summits --tempdir $(dirname $peakdir) -s 100
    


        
    
sbatch /data/public/tpadvits/PROJECTS/PodocytePJ/Chip_seq/Pipeline/XX.slurm



readdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1
samples=$(find $readdir|egrep '*Rep2_.bam$')
for i in $samples;
do
    echo $i
    name=${i/_.bam/}
    #samtools fixmate -r $i  $name"_filtered2.bam"
    samtools view --threads 7 -b -q 30 $i > $name"_MAPQ30.bam"
    samtools index  $name"_MAPQ30.bam" 
    samtools rmdup -s  $name"_MAPQ30.bam"  $name"_MAPQ30filt.bam"
    samtools index $name"_MAPQ30filt.bam"
done

samples=$(find $readdir|egrep '*_168.*_MAPQ30.bam$')
for i in $samples;
do
    name=${i/_MAPQ30.bam/}
    samtools rmdup -s  $i  $name"_MAPQ30filt.bam"
    samtools index $name"_MAPQ30filt.bam"
    samtools view -f 0x40 -b $i > $name"_MAPQ30filtSE.bam"
done

peakdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/peaks/Wt1/peaks_SEbwa_merged
bamdir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30
cd $peakdir
    
NAME="WT1_het_SE_14.06.17"
    macs2 callpeak -f BAM -t $bamdir/SN7640223_16876_ChipWt1hetF0513Rep3_MAPQ30filtSE.bam $bamdir/SN7640285_29526_ChipWt1hetI2013Rep1_MAPQ30filtSE.bam -c $bamdir/SN7640223_16877_InputWt1hetRep3_MAPQ30filtSE.bam $bamdir/SN7640285_29527_InputWt1hetRep1_MAPQ30filtSE.bam -n $NAME -B -g mm -q 0.05 --call-summits --tempdir $(dirname $peakdir)

    
## eRNA analysis
# get bed file with all intrones excluding intersections with exones and 3' 5' UTRs
dir=/data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38
cd /data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis
cat /data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis/WT1peaks_hetwt_cons.bed | bedtools sort |  bedtools subtract -a stdin -b $dir/GRCm38_GENCODE_VM11_comrehensive_exones.bed | bedtools subtract -a stdin -b $dir/GRCm38_GENCODE_VM11_comprehensive_3prim.bed | bedtools subtract -a stdin -b $dir/GRCm38_GENCODE_VM11_comprehensive_5prim.bed | bedtools subtract -a stdin -b $dir/GRCm38_GENCODE_VM11_comrehensive_intrones.bed | bedtools sort > GRCm38_WT1peaks_hetwtcons.bed


kallisto index --make-unique -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.WT1hetwtcons.idx /data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis/WT1peaks_hetwtcons_sequence.fa

kallisto index --make-unique -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.WT1hetwtconsNOcod.idx /data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis/WT1peaks_sequenceNOcoding.fa


data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/WT1/glomeruli_Ribo_minus/trimmed
index=/data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.WT1hetwtconsNOcod.idx
cd $data_dir
samples=$(find * -type f -name "*_1P")
outdir=/data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis/kallisto_WT1NOcod
for fq in $samples;
do
  fq2=${fq/_1P/_2P}
  dirname=${fq/filtered_1P/counted_dbWT1NOcod}
  echo $dirname
  kallisto quant -i $index --rf-stranded -o $outdir/$dirname -b 100 --bias -t 7 $fq $fq2 
done




kallisto index --make-unique -i /data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.cDNAncRNA_WT1.idx /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.cdna_ncRNA_WT1.fa

data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/WT1/glomeruli_Ribo_minus/trimmed
index=/data/user/tpadvits/BIOSOFT/kallisto_linux-v0.43.0/Mmus.GRCm38.cDNAncRNA_WT1.idx
cd $data_dir
samples=$(find * -type f -name "*_1P")
outdir=/data/user/tpadvits/PROJECTS/MARTIN_PJ/Chipseq_analysis/WT1/Chipseq_RNAseq_integ/eRNA_analysis/kallisto_cDNAncRNA_WT1hetwtcons
for fq in $samples;
do
  fq2=${fq/_1P/_2P}
  dirname=${fq/filtered_1P/counted_cDNAncRNA_WT1hetwt}
  echo $dirname
  kallisto quant -i $index --rf-stranded -o $outdir/$dirname -b 100 --bias -t 7 $fq $fq2 
done


##### subtract one BAM file from another #########
dir=/cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30
cd /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/Data/bam/Wt1/filtered_SE_MAPQ30/merged
bamCompare -b1 $dir/SN7640223_16876_ChipWt1hetF0513Rep3_MAPQ30filtSE.bam  -b2 $dir/SN7640223_16877_InputWt1hetRep3_MAPQ30filtSE.bam -o hetRep3_normRPKM.bw --ratio subtract --normalizeUsingRPKM
bamCompare -b1 $dir/SN7640285_29526_ChipWt1hetI2013Rep1_MAPQ30filtSE.bam  -b2 $dir/SN7640285_29527_InputWt1hetRep1_MAPQ30filtSE.bam -o hetRep1_normRPKM.bw --ratio subtract --normalizeUsingRPKM


bamCompare -b1 $dir/SN7640223_16878_ChipWt1wtF0513Rep3_MAPQ30filtSE.bam  -b2 $dir/SN7640223_16879_InputWt1wtRep3_MAPQ30filtSE.bam -o wtRep3_normRPKM.bw --ratio subtract --normalizeUsingRPKM
bamCompare -b1 $dir/SN7640223_16880_ChipWt1wtI2013Rep1_MAPQ30filtSE.bam  -b2 $dir/SN7640223_16881_InputWt1wtRep1_MAPQ30filtSE.bam -o wtRep1_normRPKM.bw --ratio subtract --normalizeUsingRPKM


### merge close features 
bedtools merge [OPTIONS] -i <BED/GFF/VCF/BAM>

### 
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703428
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703427
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703426
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703425
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703424
fastq-dump  --outdir /cellfile/DATASETS/Mouse/Kidney/Chip-seq/Bulk/In-house/ChIPseq_Wt1_Tead1_Wt1hetdel/analysis/Wt1/martinWT1_fasta -I -W  SRR1703423
