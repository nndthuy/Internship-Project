import concurrent.futures
import threading
import argparse
import itertools
import os
import gzip
import logging
from functools import partial
from dataclasses import dataclass
from typing import List, Tuple, Set

logging.basicConfig(level=logging.INFO)

global_lock = threading.Lock()

@dataclass
class Seq:
    id: str
    sequence: str

def deduplicate_chunk(sequences: List[Seq], global_results: List[Seq], seen_sequences: Set[str]) -> List[Seq]:
    local_filtered_seqs = []

    for seq_obj in sequences:
        contained = any(seq_obj.sequence in other_seq_obj.sequence for other_seq_obj in sequences if len(seq_obj.sequence) < len(other_seq_obj.sequence))
        if not contained:
            with global_lock:
                if seq_obj.sequence not in seen_sequences:
                    seen_sequences.add(seq_obj.sequence)
                    global_results.append(seq_obj)
                    local_filtered_seqs.append(seq_obj)
    return local_filtered_seqs

def recursive_deduplication(sequences: List[Seq], num_threads: int) -> List[Seq]:
    chunk_size = len(sequences) // num_threads
    seq_chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]

    global_results = []
    seen_sequences = set()

    # Using partial to freeze some arguments of deduplicate_chunk
    worker = partial(deduplicate_chunk, global_results=global_results, seen_sequences=seen_sequences)

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(worker, seq_chunks))

    combined = global_results
    before_deduplication = len(combined)
    deduplicated = deduplicate_chunk(combined, [], set())
    after_deduplication = len(deduplicated)

    # If the deduplicated result is shorter than the combined result, recurse
    if before_deduplication > after_deduplication:
        return recursive_deduplication(deduplicated, num_threads)
    else:
        return deduplicated

def read_file(file: str) -> List[str]:
    if file.endswith('.gz'):
        with gzip.open(file, 'rt') as gz_file:
            return gz_file.readlines()
    else:  # Plain text
        with open(file, 'r') as txt_file:
            return txt_file.readlines()

def process_file(input_path: str) -> List[Seq]:
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

def deduplicate_fasta(input_file: str, output_file: str, num_threads: int) -> None:
    sequences = process_file(input_file)

    sequences.sort(key=lambda seq_obj: len(seq_obj.sequence), reverse=True)

    result = recursive_deduplication(sequences, num_threads)

    with open(output_file, 'w') as outfile:
        for seq_obj in result:
            outfile.write(">" + seq_obj.id + "\n")
            outfile.write(seq_obj.sequence + "\n")

def main():
    parser = argparse.ArgumentParser(description='Deduplicate FASTA sequences.')
    parser.add_argument('-i', '--all_nematoda_fasta', help='/home/diemthuyy/Documents/INTERNSHIP/nematoda/download_seq')
    parser.add_argument('-o', '--uniq_all_nematoda_fasta', help='/home/diemthuyy/Documents/INTERNSHIP/nematoda/download_seq')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    max_threads = os.cpu_count()
    if not isinstance(args.num_threads, int) or args.num_threads < 1:
        logging.error("The number of threads must be a positive integer.")
        raise ValueError("The number of threads must be a positive integer.")
    if args.num_threads > max_threads:
        print(f"Warning: You've specified more threads ({args.num_threads}) than available CPU cores ({max_threads}). Using {max_threads} threads instead.")
        args.num_threads = max_threads

    deduplicate_fasta(args.all_nematoda_fasta, args.uniq_all_nematoda_fasta, args.num_threads)

if __name__ == '__main__':
    main()

