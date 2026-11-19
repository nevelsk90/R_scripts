## Genrich https://github.com/jsh58/Genrich 
cd /home/tim_nevelsk/SOFT/Genrich-master


#### Wt1 het.del ###
inputdir=/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Wt1.het.del/mapped.BWA
outputdir=/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Wt1.het.del/Peaks/Genrich

# samtools
ffiles=$(find /media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Wt1.het.del/mapped.BWA/PairedEnds -type f -name "*.bam")
for fq in $ffiles
    do
    samtools sort $fq -o ${fq/MAPQ30filt.sorted/MAPQ30filt} -m 16000M
#     samtools sort -n $fq -o ${fq/MAPQ30filtSE/MAPQ30filtSE.sorted} -m 16000M     # for Genrich sorted by n
    samtools index -b ${fq/MAPQ30filt.sorted/MAPQ30filt} -o ${fq/MAPQ30filt.sorted.bam/MAPQ30filt.bam.bai} --threads 4
    done

### Chipseq, mixture of SE and PE. running only on SE gives to many FP
# hetdel
./Genrich  -t "${inputdir}/SE/SN7640285_29526_ChipWt1hetI2013Rep1_MAPQ30filtSE.sorted.bam","${inputdir}/SE/SN7640285_29522_ChIPWt1hetF0513Rep2_MAPQ30filtSE.sorted.bam","${inputdir}/PairedEnds/SN7640223_16876_ChipWt1hetF0513Rep3_MAPQ30filt.sorted.bam"\
 -c "${inputdir}/SE/SN7640285_29527_InputWt1hetRep1_MAPQ30filtSE.sorted.bam","${inputdir}/SE/SN7640285_29523_InputWt1hetRep2_MAPQ30filtSE.sorted.bam","${inputdir}/PairedEnds/SN7640223_16877_InputWt1hetRep3_MAPQ30filt.sorted.bam"\
 -o "${outputdir}/Wt1hetdel.narrowPeak" -r -y -v -q 1

 # control
 ./Genrich  -t "${inputdir}/PairedEnds/SN7640223_16880_ChipWt1wtI2013Rep1_MAPQ30filt.sorted.bam","${inputdir}/SE/SN7640285_29524_ChIPWt1wtF0513Rep2_MAPQ30filtSE.sorted.bam","${inputdir}/PairedEnds/SN7640223_16878_ChipWt1wtF0513Rep3_MAPQ30filt.sorted.bam"\
 -c "${inputdir}/PairedEnds/SN7640223_16881_InputWt1wtRep1_MAPQ30filt.sorted.bam","${inputdir}/SE/SN7640285_29525_InputWt1wtRep2_MAPQ30filtSE.sorted.bam","${inputdir}/PairedEnds/SN7640223_16879_InputWt1wtRep3_MAPQ30filt.sorted.bam"\
 -o "${outputdir}/Wt1wt.narrowPeak" -r -y -v -q 1
 
 
 ### Tead1 het.del ###
inputdir=/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Tead1_Tead1.Wt1hd/mapped.BWA
outputdir=/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Tead1_Tead1.Wt1hd/Peaks/Genrich

# samtools
ffiles=$(find $inputdir -type f -name "*.bam")
for fq in $ffiles
    do
#     samtools sort $fq -o ${fq/MAPQ30filt.sorted/MAPQ30filt} -m 16000M
    samtools sort -n $fq -o ${fq/MAPQ30unique/MAPQ30unique.sorted} -m 16000M     # for Genrich sorted by n
#     samtools index -b ${fq/MAPQ30filt.sorted/MAPQ30filt} -o ${fq/MAPQ30filt.sorted.bam/MAPQ30filt.bam.bai} --threads 4
    done

### Chipseq, mixture of SE and PE. running only on SE gives to many FP
# hetdel
./Genrich  -t "${inputdir}/K002000310_108865_S17_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108873_S21_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108875_S22_L004_R1_001_MAPQ30unique.sorted.bam"\
 -c "${inputdir}/K002000310_108853_S11_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108861_S15_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108863_S16_L004_R1_001_MAPQ30unique.sorted.bam"\
 -o "${outputdir}/Tead1.Wt1hetdel.narrowPeak" -r -y -v -q 0.1

 # control
 ./Genrich  -t "${inputdir}/K002000310_108867_S18_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108869_S19_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108871_S20_L004_R1_001_MAPQ30unique.sorted.bam"\
 -c "${inputdir}/K002000310_108855_S12_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108857_S13_L004_R1_001_MAPQ30unique.sorted.bam","${inputdir}/K002000310_108859_S14_L004_R1_001_MAPQ30unique.sorted.bam"\
 -o "${outputdir}/Tead1.Wt1wt.narrowPeak" -r -y -v -q 0.1
 


 
#  ### MSPC
#  inputPdir=/media/tim_nevelsk/WD_tim/PROJECTS/PODOCYTES/Chipseq/Wt1.het.del/Peaks
#  
#  ./mspc -i ${inputPdir}/WT1.2018.peaks/WT1_het_rep1p0.01_04.01.18_peaks.narrowPeak ${inputPdir}/WT1.2018.peaks/WT1_het_rep2p0.01_04.01.18_peaks.narrowPeak ${inputPdir}/WT1.2018.peaks/WT1_het_rep3p0.01_04.01.18_peaks.narrowPeak -r bio -w 1e-5 -s 1e-10 -c 2
# 
#   ./mspc -i ${inputPdir}/WT1.2018.peaks/WT1_wt_rep1p0.01_04.01.18_peaks.narrowPeak ${inputPdir}/WT1.2018.peaks/WT1_wt_rep2p0.01_04.01.18_peaks.narrowPeak ${inputPdir}/WT1.2018.peaks/WT1_wt_rep3p0.01_04.01.18_peaks.narrowPeak -r bio -w 1e-5 -s 1e-10 -c 2
