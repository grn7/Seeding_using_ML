# script to label k-mers with their success rates (no. of times they appear in the hit window/ total no. of times they appear in reference genome)

import pickle
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from build_kmer_index import get_reverse_complement

total_hits = defaultdict(int) # dictionary to count total hits for each k-mer
total_occurrences = defaultdict(int) # dictionary to count total occurrences of each k-mer in the reference genome
reads_seen = defaultdict(set) # keeps a set of unique read IDs where the k-mer was found

# load pre-built index
print("Loading k-mer index...")
with open('k_mer_index_k_11.pkl', 'rb')as f:
    kmer_index = pickle.load(f)

def process_reads(fastq_file, truth_tsv, window_multiplier):
    print(f"Processing {fastq_file}...")
    
    # load truth table 
    truth_df = pd.read_csv(truth_tsv, sep='\t')
    # convert to dictionary for faster lookup, using read id as key and other fields as values
    truth_dict = truth_df.set_index('read_id').to_dict('index')

    for read in SeqIO.parse(fastq_file, 'fastq'): # pulls one read at a time
        read_id = read.id

        # check if we have ground truth for this read
        if read_id in truth_dict:
            truth = truth_dict[read_id] # storing all values for this read
            true_origin = truth['ref_start'] # true origin of read in reference genome
            read_length = truth['read_length'] 

            # calc hit window
            hit_window = window_multiplier * read_length
            lower_bound = max(0, true_origin - int(hit_window))
            upper_bound = true_origin + read_length + int(hit_window)
            
            read_seq = str(read.seq)
            k = 11
            for i in range(len(read_seq) - k + 1):
                kmer = read_seq[i:i+k]

                if 'N' in kmer:
                    continue

                reverse_comp_kmer = get_reverse_complement(kmer)

                if kmer <= reverse_comp_kmer:
                    canonical_kmer = kmer
                else:
                    canonical_kmer = reverse_comp_kmer

                
                positions = kmer_index.get(canonical_kmer, []) # returns list of position where k-mer is in reference genome
                # list is empty if k-mer not found in reference genome 

                if len(positions) == 0: # this occurs if k-mer is some mutation that is not present in reference genome
                    # skip entirely — error k-mer, not in reference, never in LUT
                    continue 

                # only reach here if k-mer is a real reference k-mer
                reads_seen[canonical_kmer].add(read_id) # add read ID to set of reads where this k-mer is seen

                for p in positions:
                    total_occurrences[canonical_kmer] += 1 # count total occurrences of this k-mer in reference genome
                    if lower_bound <= p <= upper_bound:
                        total_hits[canonical_kmer] += 1 # count hits for this k-mer if it falls within the hit window


# process short reads 
process_reads('02_simulated_short.fastq', '02_short_read_truth.tsv', 1.0)

# process long reads
process_reads('02_simulated_long_0001.fastq', '02_long_read_truth.tsv', 0.5)

print("Filtering and writing final dataset\n")
# write everything to a single combined file
with open('labeled_dataset_k11.csv', 'w') as f:
    f.write('canonical_kmer,success_rate,observation_count\n')
    for kmer in reads_seen.keys():
        obs_count = len(reads_seen[kmer])
        if obs_count >= 5:
            success_rate = total_hits[kmer] / total_occurrences[kmer]
            f.write(f"{kmer},{success_rate:.4f},{obs_count}\n")

print("Dataset Created")