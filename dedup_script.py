
import concurrent.futures
import argparse
import os
import gzip
import logging
import hashlib
from functools import partial
from dataclasses import dataclass
from typing import List, Set

logging.basicConfig(level=logging.INFO)

@dataclass
class Seq:
    id: str
    sequence: str

def read_file(file: str) -> List[str]:
    """Read a file. Supports gzipped files."""
    if file.endswith('.gz'):
        with gzip.open(file, 'rt') as gz_file:
            return gz_file.readlines()
    else:
        with open(file, 'r') as txt_file:
            return txt_file.readlines()

def process_file(input_path: str) -> List[Seq]:
    """Process a FASTA file and return a list of Seq objects."""
    try:
        file_lines = read_file(input_path)
    except Exception as e:
        logging.error(f"An error occurred while reading the file: {str(e)}")
        return []

    seqs_with_id = []
    cur_id = None
    cur_seq = ""
    for line in file_lines:
        line = line.strip()
        if line.startswith(">"):
            if cur_id and cur_seq:
                seqs_with_id.append(Seq(cur_id, cur_seq))
            cur_id = line[1:]
            cur_seq = ""
        else:
            cur_seq += line
    if cur_id and cur_seq:
        seqs_with_id.append(Seq(cur_id, cur_seq))
    return seqs_with_id

def hash_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()

def deduplicate_chunk_combined(seq_chunk: List[Seq], prior_sequences: List[Seq], seen_hashes: Set[str]) -> List[Seq]:
    """Deduplicate sequences using both hashing and containment checks."""
    unique_seqs = []
    for current_seq in seq_chunk:
        seq_hash = hash_sequence(current_seq.sequence)

        # Exact match deduplication
        if seq_hash in seen_hashes:
            continue

        is_contained = any(current_seq.sequence in prior_seq.sequence for prior_seq in prior_sequences)
        if not is_contained:
            unique_seqs.append(current_seq)
            seen_hashes.add(seq_hash)
    return unique_seqs

def recursive_deduplication_combined(sequences: List[Seq], num_threads: int) -> List[Seq]:
    """Recursively deduplicate sequences using multiple threads and combined technique."""
    chunk_size = len(sequences) // num_threads
    seq_chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]

    deduped_results = []
    prior_seqs = []
    seen_hashes = set()

    deduplication_function = partial(deduplicate_chunk_combined, prior_sequences=prior_seqs, seen_hashes=seen_hashes)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        for seq_chunk in seq_chunks:
            deduped_chunk = executor.submit(deduplication_function, seq_chunk).result()
            deduped_results.extend(deduped_chunk)
            prior_seqs.extend(seq_chunk)

    return deduped_results

def deduplicate_fasta(input_file: str, output_file: str, num_threads: int) -> None:
    """Main deduplication function using the combined technique."""
    sequences = process_file(input_file)
    sequences.sort(key=lambda seq_obj: len(seq_obj.sequence), reverse=True)

    result = recursive_deduplication_combined(sequences, num_threads)

    with open(output_file, 'w') as outfile:
        for seq_obj in result:
            outfile.write(">" + seq_obj.id + "\n")
            outfile.write(seq_obj.sequence + "\n")

def main():
    parser = argparse.ArgumentParser(description='Deduplicate FASTA sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='/home/diemthuyy/Documents/INTERNSHIP/nematoda/download_seq')
    parser.add_argument('-o', '--output_file', required=True, help='/home/diemthuyy/Documents/INTERNSHIP/nematoda/download_seq')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    deduplicate_fasta(args.input_file, args.output_file, args.num_threads)

if __name__ == '__main__':
    main()
