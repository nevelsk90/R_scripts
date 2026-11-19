#!/bin/bash -l

#SBATCH --account=tpadvits
#SBATCH --job-name=nfcore.rnaseq_glom.circ
#SBATCH --ntasks=12
#SBATCH --mem=50G
#SBATCH --error=/cellfile/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/nf-core_RNAseq/nfcore.RNAseq.error
#SBATCH --output=/cellfile/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/nf-core_RNAseq/nfcore.RNAseq.out

conda activate /cellfile/ypaul1/my-envs/nfcore_rnaseq

nextflow run nf-core/rnaseq \
	-r 3.18.0 \
	-profile docker \
	--input "/cellfile/DATASETS/Mouse/Kidney/Transcriptomics/Bulk/In_house/bulkRNA_glomCircadian_KatharinaLemberg/Data/fastq/input_rnaseq.csv" \
	--outdir "/cellfile/tpadvits/PROJECTS/PodocytePJ/RNA_seq/Public_data/Mouse_data/glomerular_circadian_PRJEB67585/nf-core_RNAseq/" \
  	--gtf "/cellfile/ypaul1/genome_assemblies/Mus_musculus.GRCm39.113.gtf.gz" \
	--fasta "/cellfile/ypaul1/genome_assemblies/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
	--skip_trimming true \
	--skip_linting true \
	--skip_qc true \
        -resume
        
###	--star_index  "" \
###	--salmon_index  "" \
### 	--contaminant_screening 'kraken2_bracken' \
###  	--kraken_db "/cellfile/tpadvits/global_data/k2_pluspfp_16gb_20241228" \
        
### nextflow clean shrivelled_avogadro -f
