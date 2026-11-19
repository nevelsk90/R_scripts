#!/bin/bash -l
#SBATCH --account=tpadvits
#SBATCH --cpus-per-task=8
#SBATCH --error=/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Bulk/In_house/bulkRNA_wtPodo_MahdiehRahmatollahi/Data/processed/misc/STARbulkEx.error
#SBATCH --output=/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Bulk/In_house/bulkRNA_wtPodo_MahdiehRahmatollahi/Data/processed/misc/STARbulkEx.out
#SBATCH --workdir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/wtPodocytes/
#SBATCH --job-name=STARbulkInt
#SBATCH --partition=all
#SBATCH --ntasks=1


#### 
# module load STAR-2.6.0c
# module load samtools-1.4.0
# module unload python3 #necessary for correct work of MACS2 and other soft
# module load python #necessary for correct work of MACS2 and other soft

#module unload bwa-0.7.7-haswell
# module load MACS2-2.1.1-20160309 #14.05.18 WT1martin
#module load MACS2-2.0.10.20131216

#module load R-3.6.0

#CMD="/cellfile/datapublic/tpadvits/PROJECTS/utility.SCRIPTS/Pipelines/bwaANDqc_raw.sh"
 bash /cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/DATA/RNA_seq/STAR2.7bulk_mapp.sh
# bash /cellfile/datapublic/tpadvits/PROJECTS/GRNinference/podocyte_net/MiraldiWay/priorGenerator_MiraldiGRN.sh
# bash /cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/DATA/RNA_seq/QC_kallistoBUS_scRNAseq.sh
# bash /cellfile/datapublic/tpadvits/PROJECTS/GRNinference/podocyte_net/MiraldiWay/Genrich.sh
# setenv("R_LIBS_USER","/cellfile/datapublic/tpadvits/SOFT/RLibs/R-3.4.2/")
# module load R
# module use /data/public/tpadvits/global_data/MyModules
# module load R/R-3.4.2tpadvits
# which R

# module load MACS2-2.0.10.20131216
# module unload python
# module unload python3
# 
# module load python

#module load python3


#module load STAR-2.5.2b
#module load python
#module load HTSeq-0.6.1p1
#module load bowtie2-2.2.9
#module load bwa-0.7.7-haswell 

#bash $CMD
#Rscript /cellfile/DATASETS/Mouse/Kidney/Chip_seq/Bulk/In_house/ChIPseq_Tead1/analysis/genoGAM/gGAM_pvalTOpeaks.R
