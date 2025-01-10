import os
import glob
from multiprocessing import Pool, Manager

# Define the input and output file paths
sam_file = '../../../../../jmolinav/Eliminar_reads_solapantes_A12/Results/solap_A12_reg_L002_ALIGNED_filtered.sam'  # The SAM file
fastq_dir = '../2/'  # Directory containing the FASTQ files
output_dir = './'  # Directory for filtered output files

# Step 1: Read identifiers from the SAM file in streaming mode
def get_identifiers(sam_file):
    print(f"[INFO] Reading identifiers from the SAM file: {sam_file}")
    identifiers = set()
    line_count = 0
    try:
        with open(sam_file, 'r') as sam:
            for line in sam:
                line_count += 1
                if line_count % 1000000 == 0:
                    print(f"[DEBUG] Processed {line_count} lines from the SAM file")
                if not line.startswith('@'):  # Ignore headers
                    fields = line.split('\t')
                    identifiers.add(fields[0])  # Identifier is in the first column
                # Additional diagnostics to check the growth of the set
                if line_count % 1000000 == 0:
                    print(f"[DEBUG] Current size of the identifier set: {len(identifiers)}")
    except MemoryError:
        print(f"[ERROR] Ran out of memory while processing line {line_count}")
        raise
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred on line {line_count}: {e}")
        raise
    print(f"[INFO] Finished reading the SAM file. Total lines processed: {line_count}. Unique identifiers: {len(identifiers)}.")
    return identifiers

# Step 2: Process a single FASTQ file in streaming mode
def process_fastq(args):
    fastq_file, identifiers = args
    output_fastq = os.path.join(output_dir, os.path.basename(fastq_file).replace('.fastq', '_filtered.fastq'))
    print(f"[INFO] Processing FASTQ file: {fastq_file}")
    try:
        with open(fastq_file, 'r') as fastq, open(output_fastq, 'w') as output:
            while True:
                header = fastq.readline().strip()  # Line 1: header
                if not header:  # End of file
                    break
                sequence = fastq.readline().strip()  # Line 2: sequence
                plus = fastq.readline().strip()     # Line 3: '+' symbol
                quality = fastq.readline().strip()  # Line 4: quality

                fastq_id = header.split()[0][1:]  # Extract identifier
                if fastq_id not in identifiers:  # Filter by identifiers
                    output.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
    except Exception as e:
        print(f"[ERROR] Error while processing FASTQ file {fastq_file}: {e}")
        raise
    print(f"[INFO] Filtering completed for: {fastq_file} -> {output_fastq}")
    return f"Filtering completed for: {fastq_file} -> {output_fastq}"

# Step 3: Prepare and parallelize
if __name__ == '__main__':
    print("[INFO] Starting the script...")
    try:
        # Get identifiers from the SAM file
        identifiers = get_identifiers(sam_file)

        # Check if the output directory exists
        if not os.path.exists(output_dir):
            print(f"[INFO] Creating the output directory: {output_dir}")
            os.makedirs(output_dir)

        # Get a list of FASTQ files
        fastq_files = glob.glob(os.path.join(fastq_dir, '*.fastq'))
        print(f"[INFO] Found {len(fastq_files)} FASTQ files in {fastq_dir}.")

        # Use a manager to share identifiers between processes
        with Manager() as manager:
            shared_identifiers = manager.list(identifiers)  # Share the set as a list
            with Pool(processes=4) as pool:
                print("[INFO] Starting parallel processing...")
                results = pool.map(process_fastq, [(f, shared_identifiers) for f in fastq_files])

            # Display processing results
            print("[INFO] Processing completed. Results:")
            for result in results:
                print(result)
    except MemoryError:
        print("[ERROR] The script ran out of memory during execution.")
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred: {e}")
