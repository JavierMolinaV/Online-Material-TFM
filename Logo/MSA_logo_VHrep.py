from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Leer las secuencias desde un archivo de texto
with open("MSA_LOGO_VHrep/longCDR3_JH1_aa1_sequences.txt", "r") as file:
    sequences = file.readlines()

# Eliminar los saltos de línea (y cualquier espacio extra)
sequences = [seq.strip() for seq in sequences]

# Convertir las secuencias a objetos SeqRecord
seq_records = [SeqRecord(Seq(seq), id=f"seq{i}") for i, seq in enumerate(sequences)]

# Guardar las secuencias en un archivo FASTA temporal
with open("temp_sequences.fasta", "w") as fasta_file:
    from Bio import SeqIO
    SeqIO.write(seq_records, fasta_file, "fasta")

# Usar ClustalW para hacer el alineamiento
clustalw_exe = "/usr/bin/clustalw"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="temp_sequences.fasta", gapopen=10, gapext=1)
stdout, stderr = clustalw_cline()

# Leer el archivo de salida del alineamiento
alignment = AlignIO.read("temp_sequences.aln", "clustal")
aligned_sequences = [str(record.seq) for record in alignment]

# Calcular la frecuencia de cada aminoácido en cada columna
all_chars = list("ACDEFGHIKLMNPQRSTVWY-")  # Incluye los 20 aminoácidos y el gap '-'
seq_matrix = np.array([list(seq) for seq in aligned_sequences])

freq_matrix = []
for i in range(seq_matrix.shape[1]):  # Iterar por cada columna
    column = seq_matrix[:, i]
    freq = {char: np.sum(column == char) / len(column) for char in all_chars}
    freq_matrix.append([freq[char] for char in all_chars])

# Convertir a un DataFrame para usarlo en logomaker
freq_matrix_df = pd.DataFrame(freq_matrix, columns=all_chars)

# Usar logomaker para generar el logo
logo = logomaker.Logo(freq_matrix_df)

# Filtrar columnas de gaps (eliminamos '-')
all_chars_no_gaps = list("ACDEFGHIKLMNPQRSTVWY")  # Excluye '-'

# Calcular frecuencia de cada aminoácido sin incluir gaps
freq_matrix = []
for i in range(seq_matrix.shape[1]):  # Iterar por cada columna
    column = seq_matrix[:, i]
    freq = {char: np.sum(column == char) / len(column) for char in all_chars_no_gaps}
    freq_matrix.append([freq[char] for char in all_chars_no_gaps])

# Convertir a un DataFrame de pandas sin gaps
freq_matrix_df = pd.DataFrame(freq_matrix, columns=all_chars_no_gaps)

# Colores personalizados para los aminoácidos
color_dict = {
    'A': '#F5CBA7', 'V': '#F5CBA7', 'L': '#F5CBA7', 'I': '#F5CBA7', 'M': '#F5CBA7', 'F': '#F5CBA7',
    'W': '#F5CBA7', 'G': '#F5CBA7', 'P': '#F5CBA7',  # Apolares: beige y ocre

    'D': '#E74C3C', 'E': '#E74C3C',                # Ácidos: rojo anaranjado

    'K': '#85C1E9', 'R': '#85C1E9', 'H': '#85C1E9', # Básicos: azul pastel y azul medio

    'N': '#A9DFBF', 'Q': '#A9DFBF', 'S': '#A9DFBF', 'T': '#A9DFBF', # Polares: verdes suaves
    'Y': '#A9DFBF', 'C': '#A9DFBF'                # Polares adicionales: gris claro
}

# Crear el logo sin gaps
logo = logomaker.Logo(freq_matrix_df , color_scheme=color_dict)

# Configurar título y mostrar el gráfico
logo.ax.set_title("Sequence Logo sin Gaps")
plt.show()



