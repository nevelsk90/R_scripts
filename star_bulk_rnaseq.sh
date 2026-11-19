#!/usr/bin/env bash
set -euo pipefail


## Run once for each genome/annotation combination:
#cd /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/homo_sapiens
#
#mkdir -p STAR_index/Homo_sapiens.GRCh38.107
#
# /cellfile/datapublic/tpadvits/SOFT/STAR-master/source/STAR \
#    --runMode genomeGenerate \
#    --runThreadN 12 \
#    --genomeDir STAR_index/Homo_sapiens.GRCh38.107 \
#    --genomeFastaFiles /cellfile/genome_assemblies/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#    --sjdbGTFfile Homo_sapiens.GRCh38.107.gtf \
#    --sjdbOverhang 149

# # cd /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/mus_musculus
# # 
# # # mkdir -p STAR_index/Mus_musculus.GRCm39
# # 
# /cellfile/datapublic/tpadvits/SOFT/STAR-master/source/STAR \
#    --runMode genomeGenerate \
#    --runThreadN 12 \
#    --genomeDir /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/mus_musculus/STAR_index/Mus_musculus.GRCm39 \
#    --genomeFastaFiles /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/mus_musculus/Mus_musculus.GRCm39.dna.toplevel.fa \
#    --sjdbGTFfile /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/mus_musculus/Mus_musculus.GRCm39.107.gtf \
#    --sjdbOverhang 149


# Usage:
#   bash star_bulk_rnaseq.sh samples.tsv /path/to/STAR_index output_dir [threads] [count_column]
#
# count_column:
#   2 = unstranded
#   3 = forward-stranded
#   4 = reverse-stranded
#
# Example:
#   bash star_bulk_rnaseq.sh samples.tsv STAR_index star_out 12 2

#Tab-separated samples.tsv:
#sample	r1	r2
#sample1	/path/sample1_R1.fastq.gz	/path/sample1_R2.fastq.gz
#sample2	/path/sample2.fastq.gz	NA

SAMPLES=${1:?Provide samples.tsv}
GENOME_DIR=${2:?Provide STAR index directory}
OUTDIR=${3:?Provide output directory}
THREADS=${4:-8}
COUNT_COL=${5:-2}

mkdir -p "${OUTDIR}"/{bam,counts,logs}

# Skip header.
tail -n +2 "$SAMPLES" |
  while IFS=$'\t' read -r sample r1 r2; do

echo "Processing ${sample}"

[[ -f "$r1" ]] || {
  echo "ERROR: R1 not found: $r1" >&2
  exit 1
}

reads=( "$r1" )

if [[ -n "${r2:-}" && "$r2" != "NA" && "$r2" != "." ]]; then
[[ -f "$r2" ]] || {
  echo "ERROR: R2 not found: $r2" >&2
  exit 1
}
reads+=( "$r2" )
fi

# Use zcat only for compressed FASTQ files.
read_command=()
if [[ "$r1" == *.gz ]]; then
read_command=( --readFilesCommand zcat )
fi

prefix="${OUTDIR}/${sample}_"

/cellfile/datapublic/tpadvits/SOFT/STAR-master/source/STAR \
--runThreadN "$THREADS" \
--genomeDir "$GENOME_DIR" \
--readFilesIn "${reads[@]}" \
"${read_command[@]}" \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard \
--quantMode GeneCounts \
--outFileNamePrefix "$prefix" \
--limitBAMsortRAM 60014287022

mv "${prefix}Aligned.sortedByCoord.out.bam" \
"${OUTDIR}/bam/${sample}.bam"

mv "${prefix}Log.final.out" \
"${OUTDIR}/logs/${sample}.Log.final.out"

# Remove STAR summary rows and select gene ID plus requested count column.
awk -v col="$COUNT_COL" '
        $1 !~ /^N_/ {print $1 "\t" $col}
    ' "${prefix}ReadsPerGene.out.tab" \
> "${OUTDIR}/counts/${sample}.counts.tsv"

# rm -f \
# #        "${prefix}ReadsPerGene.out.tab" \ # don't remove original count
# "${prefix}Log.out" \
# "${prefix}Log.progress.out" \
# "${prefix}SJ.out.tab"
# 
done
# 
# # Confirm that all count files contain genes in the same order.
# first_file=$(find "${OUTDIR}/counts" -name "*.counts.tsv" | sort | head -n 1)
# 
# [[ -n "$first_file" ]] || {
#   echo "ERROR: no count files generated" >&2
#   exit 1
# }
# 
# cut -f1 "$first_file" > "${OUTDIR}/counts/gene_ids.tmp"
# 
# for file in $(find "${OUTDIR}/counts" -name "*.counts.tsv" | sort); do
# cut -f1 "$file" > "${OUTDIR}/counts/current_gene_ids.tmp"
# 
# cmp -s \
# "${OUTDIR}/counts/gene_ids.tmp" \
# "${OUTDIR}/counts/current_gene_ids.tmp" || {
#   echo "ERROR: inconsistent gene order in $file" >&2
#   exit 1
# }
# done
# 
# # Build the matrix header.
# {
#   printf "gene_id"
#   for file in $(find "${OUTDIR}/counts" -name "*.counts.tsv" | sort); do
#   sample=$(basename "$file" .counts.tsv)
#   printf "\t%s" "$sample"
#   done
#   printf "\n"
# } > "${OUTDIR}/gene_counts.tsv"
# 
# # Combine gene IDs and count columns.
# count_columns=()
# for file in $(find "${OUTDIR}/counts" -name "*.counts.tsv" | sort); do
# tmp="${file}.column.tmp"
# cut -f2 "$file" > "$tmp"
# count_columns+=( "$tmp" )
# done
# 
# paste \
# "${OUTDIR}/counts/gene_ids.tmp" \
# "${count_columns[@]}" \
# >> "${OUTDIR}/gene_counts.tsv"
# 
# rm -f \
# "${OUTDIR}/counts/gene_ids.tmp" \
# "${OUTDIR}/counts/current_gene_ids.tmp" \
# "${OUTDIR}"/counts/*.column.tmp

echo "Finished"
# echo "Count matrix: ${OUTDIR}/gene_counts.tsv"
echo "BAM files:    ${OUTDIR}/bam/"
echo "STAR logs:    ${OUTDIR}/logs/"
