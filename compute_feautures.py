# this file computes additional features for each k-mer and prepares the dataset to be fed to ML model

import pandas as pd
import pickle
import math 

# function to calculate GC content
def gc_content(kmer):

    kmer = kmer.upper()
    g_count = kmer.count('G')
    c_count = kmer.count('C')

    gc_content = (g_count + c_count) / len(kmer)
    return gc_content

# function to calculate shannon entropy
def shannon_entropy(kmer):

    k = len(kmer)
    kmer = kmer.upper()

    count = [0] * 4 # A is 0, C 1, G 2, T 3
    count[0] = kmer.count('A')
    count[1] = kmer.count('C')
    count[2] = kmer.count('G')
    count[3] = kmer.count('T')

    prob = [0] * 4
    for i in range(len(count)):
        prob[i] = count[i] / k

    shannon_entropy = 0.0
    for i in range(len(count)):
        if count[i] != 0 :
            shannon_entropy += prob[i] * math.log2(prob[i])

    shannon_entropy = -shannon_entropy
    return shannon_entropy

# function to calculate homopolymer fraction (measure of longest run of consecutive identical bases in the kmer)
def homopolymer_fraction(kmer):

    k = len(kmer)
    kmer = kmer.upper()
    counter = set()
    count = 1

    for i in range(1, len(kmer)):
        if kmer[i] == kmer[i-1]:
            count += 1
        else:
            counter.add(count)
            count = 1
    counter.add(count) # to add the count for homopolymer which goes from a given element to last element

    homo_frac = max(counter) / k
    return homo_frac

# function to extract all consecutive 2-base substrings from k-mer and check how many are distinct
def dinucleotide_complexity(kmer):

    kmer = kmer.upper()
    dinucleotide_set = set()

    for i in range(len(kmer) - 1):
        dinucleotide_set.add(kmer[i] + kmer[i+1])

    dinucleotide_complexity = len(dinucleotide_set) / 16
    return dinucleotide_complexity

# load the kmer index and labelled k-mers
with open('k_mer_index_k_11.pkl', 'rb') as f:
    kmer_index = pickle.load(f)
df = pd.read_csv('labeled_dataset_k11.csv')

df['gc_content'] = df['canonical_kmer'].apply(gc_content)
df['entropy'] = df['canonical_kmer'].apply(shannon_entropy)
df['homopoly_frac'] = df['canonical_kmer'].apply(homopolymer_fraction)
df['dinuc_complexity'] = df['canonical_kmer'].apply(dinucleotide_complexity)

# calculate log frequency 
df['log_freq'] =df['canonical_kmer'].apply(lambda x: math.log2(len(kmer_index.get(x,[])) + 1))

# reorder columns
df = df[['canonical_kmer', 'log_freq', 'entropy', 'gc_content', 'homopoly_frac', 'dinuc_complexity', 'success_rate', 'observation_count']]

# save dataset
df.to_csv('training_data_k11.csv', index=False)









    