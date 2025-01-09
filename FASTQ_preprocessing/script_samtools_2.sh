#!/usr/bin/env bash

#$ -P ARR
#$ -A jmolinav
#$ -N samtools_2_A12
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

# Archivos de entrada
FILTER_SAM="../Results/solap_A12_reg_L002_ALIGNED_filtered.sam"
INPUT_SAM="../Results/solap_A12_reg_L002.sam"
OUTPUT_FASTQ="~/NetVolumes/data_arr/arodriguezr/jmolinav/Eliminar_reads_solapantes_A12/Results_Filter2"

# Paso 1: Extraer nombres de lecturas a descartar
samtools view "$FILTER_SAM" | awk '{print $1}' > discard_reads_2.txt

# Paso 2: Filtrar las lecturas de INPUT_SAM que no están en discard_reads.txt
samtools view -h "$INPUT_SAM" | grep -vFf discard_reads_2.txt | samtools view -Sb - > filtered_output_2.bam

# Paso 3: Convertir el BAM filtrado a FASTQ
samtools fastq filtered_output_2.bam -1 "$OUTPUT_FASTQ" > "$OUTPUT_FASTQ"

# Limpieza (opcional)
rm discard_reads_2.txt
