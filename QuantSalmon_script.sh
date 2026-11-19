# # run salmon from a docker
 # # buid indeces



docker run -it -v /cellfile/datapublic/tpadvits/ANNOTATIONS/Genome_annot/mus_musculus:/index -v /cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/DATA/RNA_seq/bulk_WT1het.del./glomeruli_Ribo_minus/FASTQ:/data -v /cellfile/datapublic/tpadvits/PROJECTS/PodocytePJ/DATA/RNA_seq/bulk_WT1het.del./glomeruli_Ribo_minus/Salmon:/out combinelab/salmon


chmod -R 777 /index


# #!/bin/bash
for fn in /data/*_1*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i /index/Mus_musculus.GRCm39_salmon.index -l A \
         -1 /data/$samp \
         -2 /data/${samp/_1/_2} \
         -p 16 --validateMappings -o /out/${samp/_1*/_quant}
done
