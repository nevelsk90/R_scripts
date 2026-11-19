#!/usr/bin/bash
   
PATH=$PATH:~/data/public/tpadvits/SOFT/sratoolkit.3.0.2-ubuntu64/bin

## for sc data where biological and seuqncing reads are in seperate SRR files 
prefetch --max-size 1t $1 
fasterq-dump  -v -e 20 $1
 
##  navigate to a destination folder, put SRR_Acc_List.txt there, run this in a terminal:
# cat SRR_Acc_List.txt | xargs -n 1 bash /cellfile/datapublic/tpadvits/PROJECTS/utility.SCRIPTS/get_SRR_data.sh 

