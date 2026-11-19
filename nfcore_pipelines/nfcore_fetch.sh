#!/bin/bash -l

#SBATCH --account=tpadvits
#SBATCH --job-name=fetch_glom.circ
#SBATCH --mem=50G
#SBATCH --ntasks=12
#SBATCH --error=fetch.error
#SBATCH --output=fetch.out

## conda activate /cellfile/ypaul1/my-envs/nfcore_rnaseq

nextflow run nf-core/fetchngs \
   -r 1.12.0 \
   -profile docker \
   --input "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Bulk/Public/glomerular_circadian_PRJEB67585/Data/ids.csv" \
   --outdir "/cellfile/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/FASTA" \
   -resume
