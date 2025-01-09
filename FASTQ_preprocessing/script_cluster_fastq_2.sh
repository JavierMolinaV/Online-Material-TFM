#!/usr/bin/env bash

#$ -P ARR
#$ -A jmolinav
#$ -N fastq_A12
#$ -V
#$ -pe pthreads 4
#$ -cwd
#$ -l h_vmem=30G
#$ -o /data_lab_ARR/arodriguezr/jmolinav/RNAvel_SP/outs
#$ -e /data_lab_ARR/arodriguezr/jmolinav/RNAvel_SP/error
# #$ -S "/usr/bin/env bash"

# Activar el entorno conda
source ~/miniconda3/bin/activate NGS

# Ruta a los ejecutables

BWA_PATH=~/miniconda3/envs/NGS/bin/
SAMTOOLS_PATH=~/miniconda3/envs/NGS/bin/

# Exportar las rutas para que los programas sepan dónde encontrar las dependencias
export PATH=$BWA_PATH:$SAMTOOLS_PATH:$PATH


~/miniconda3/envs/NGS/bin/python3.10 ../../../../../jmolinav/Eliminar_reads_solapantes_A12/script/script_A12_fastq_curation_2.py
