#!/usr/bin/env bash

#$ -P ARR
#$ -A jmolinav
#$ -N bwa_A12
#$ -V
#$ -pe pthreads 4
#$ -cwd
#$ -l h_vmem=5G
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


~/miniconda3/envs/NGS/bin/bwa mem -t4 ../Ref/solap_A12.fa ../Data/BM_VDJ_S2_L002_R2_001.fastq > ./solap_A12_reg_L002.sam

awk '{ 
    # Search for the CIGAR string in the 6th column of the SAM file
    if ($6 ~ /[0-9]+S/) {
        # Extract the number preceding 'S' using a regular expression
        num_s = substr($6, 1, index($6, "S") - 1);
        if (num_s <= 30) {
            print $0;  # Print the line if the number before 'S' is ≤ 30
        }
    } else {
        print $0;  # Print the line if the CIGAR string does not contain 'S'
    }
}' solap_A12_reg_L002.sam > solap_A12_reg_L002_filtered.sam
