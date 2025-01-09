#!/usr/bin/env bash

#$ -P ARR
#$ -A jmolinav
#$ -N A12_BM_VDJ_new
#$ -V
#$ -pe pthreads 20
#$ -cwd
#$ -l h_vmem=5G
# #$ -o /data_lab_ARR/arodriguezr/A12_10x/out/
# #$ -e /data_lab_ARR/arodriguezr/A12_10x/error/
# #$ -S "/usr/bin/env bash"


/data_lab_ARR/arodriguezr/2.Cell_ranger/Data/cellranger-7.1.0/bin/cellranger count --id=A12_BM --transcriptome /data_lab_ARR/arodriguezr/2.Cell_ranger/Data/New_ref_A12_16_09_24 --fastqs /data_lab_ARR/arodriguezr/2.Cell_ranger/Data_A12/BM_VDJ --sample BM_VDJ 
