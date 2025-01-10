import os
import glob
from multiprocessing import Pool, Manager
# Definir las rutas de los archivos de entrada y salida
sam_file = '../../../../../jmolinav/Eliminar_reads_solapantes_A12/Results/solap_A12_reg_L001_ALIGNED_filtered.sam'  # El archivo SAM
fastq_dir = '../1/'  # El directorio donde se encuentran los archivos FASTQ
output_dir = './'  # Directorio para los archivos filtrados

# Paso 1: Leer identificadores del archivo SAM en streaming
def get_identifiers(sam_file):
    print(f"[INFO] Leyendo identificadores del archivo SAM: {sam_file}")
    identifiers = set()
    line_count = 0
    try:
        with open(sam_file, 'r') as sam:
            for line in sam:
                line_count += 1
                if line_count % 1000000 == 0:
                    print(f"[DEBUG] Procesadas {line_count} líneas del archivo SAM")
                if not line.startswith('@'):  # Ignorar encabezados
                    fields = line.split('\t')
                    identifiers.add(fields[0])  # Identificador en la primera columna
                # Diagnóstico adicional para verificar el crecimiento del conjunto
                if line_count % 1000000 == 0:
                    print(f"[DEBUG] Tamaño actual del conjunto de identificadores: {len(identifiers)}")
    except MemoryError:
        print(f"[ERROR] Se quedó sin memoria al procesar la línea {line_count}")
        raise
    except Exception as e:
        print(f"[ERROR] Ocurrió un error inesperado en la línea {line_count}: {e}")
        raise
    print(f"[INFO] Se completó la lectura del archivo SAM. Total de líneas procesadas: {line_count}. Identificadores únicos: {len(identifiers)}.")
    return identifiers

# Paso 2: Procesar un único archivo FASTQ en streaming
def process_fastq(args):
    fastq_file, identifiers = args
    output_fastq = os.path.join(output_dir, os.path.basename(fastq_file).replace('.fastq', '_filtered.fastq'))
    print(f"[INFO] Procesando archivo FASTQ: {fastq_file}")
    try:
        with open(fastq_file, 'r') as fastq, open(output_fastq, 'w') as output:
            while True:
                header = fastq.readline().strip()  # Línea 1: encabezado
                if not header:  # Fin del archivo
                    break
                sequence = fastq.readline().strip()  # Línea 2: secuencia
                plus = fastq.readline().strip()     # Línea 3: símbolo '+'
                quality = fastq.readline().strip()  # Línea 4: calidad

                fastq_id = header.split()[0][1:]  # Extraer identificador
                if fastq_id not in identifiers:  # Filtrar por identificadores
                    output.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    except Exception as e:
        print(f"[ERROR] Error al procesar el archivo FASTQ {fastq_file}: {e}")
        raise
    print(f"[INFO] Filtrado completado para: {fastq_file} -> {output_fastq}")
    return f"Filtrado completado para: {fastq_file} -> {output_fastq}"

# Paso 3: Preparar y paralelizar
if __name__ == '__main__':
    print("[INFO] Iniciando el script...")
    try:
        # Obtener identificadores del archivo SAM
        identifiers = get_identifiers(sam_file)

        # Verificar si el directorio de salida existe
        if not os.path.exists(output_dir):
            print(f"[INFO] Creando el directorio de salida: {output_dir}")
            os.makedirs(output_dir)

        # Obtener lista de archivos FASTQ
        fastq_files = glob.glob(os.path.join(fastq_dir, '*.fastq'))
        print(f"[INFO] Se encontraron {len(fastq_files)} archivos FASTQ en {fastq_dir}.")

        # Usar un gestor para compartir identificadores entre procesos
        with Manager() as manager:
            shared_identifiers = manager.list(identifiers)  # Compartir el conjunto como lista
            with Pool(processes=4) as pool:
                print("[INFO] Iniciando procesamiento en paralelo...")
                results = pool.map(process_fastq, [(f, shared_identifiers) for f in fastq_files])

            # Mostrar resultados del procesamiento
            print("[INFO] Procesamiento completado. Resultados:")
            for result in results:
                print(result)
    except MemoryError:
        print("[ERROR] El script se quedó sin memoria durante la ejecución.")
    except Exception as e:
        print(f"[ERROR] Ocurrió un error inesperado: {e}")
