#!/usr/bin/env bash

#$ -P ARR
#$ -A arodriguezr
#$ -N A12_exp_10x_SP_ClusterJob
#$ -V
#$ -pe pthreads 20
#$ -cwd
#$ -l h_vmem=5G
# #$ -o /data_lab_ARR/arodriguezr/A12_10x/out/
# #$ -e /data_lab_ARR/arodriguezr/A12_10x/error/
# #$ -S "/usr/bin/env bash"


/data_lab_ARR/arodriguezr/2.Cell_ranger/Data/cellranger-7.1.0/bin/cellranger multi --id=A12_SP --csv=/data_lab_ARR/arodriguezr/2.Cell_ranger/Data_A12/multi_config_A12_exp_10x_SP_mod_ref.csv --jobmode=local --localcores=4
