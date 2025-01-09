#!/usr/bin/env bash

#$ -P ARR
#$ -A jmolinav
#$ -N RNAvel
#$ -V
#$ -pe pthreads 20
#$ -cwd
#$ -l h_vmem=5G
#$ -o /data_lab_ARR/arodriguezr/jmolinav/RNAvel/outs
#$ -e /data_lab_ARR/arodriguezr/jmolinav/RNAvel/error
# #$ -S "/usr/bin/env bash"

# Activar el entorno conda
source ~/miniconda3/bin/activate RNAvel

# Ruta a los ejecutables

VELOCYTO_PATH=~/miniconda3/envs/RNAvel/bin/
SAMTOOLS_PATH=~/miniconda3/envs/RNAvel/bin/

# Exportar las rutas para que los programas sepan dónde encontrar las dependencias
export PATH=$VELOCYTO_PATH:$SAMTOOLS_PATH:$PATH


~/miniconda3/envs/RNAvel/bin/velocyto run -b ./barcodes_BM.tsv -o Results_BM -m ./mm10_rmsk.gtf ./sample_alignments_BM.bam ./Mus_musculus.GRCm39.112.chr.gtf
