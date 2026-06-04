# python code to build k-mer index from a fasta file

from Bio import SeqIO # SeqIO - Sequence Input/Output. Designed to write and read biological data formats fast
from collections import defaultdict # used to immdiately append position of a k-mer to a key without checking if the key exists
import pickle # to save k-mer index to a file for later use 

# function to get reverse complement of a DNA sequence
def get_reverse_complement(seq): # take a string(k-mer) as input

    # define the complement using a dictionary
    complement = {'A': 'T',
                  'T': 'A',
                  'C': 'G',
                  'G': 'C'}
    
    # intialize empty string to store reverse complement
    reverse_complement = []
    for base in seq: # base is a char
        comp_base = complement.get(base.upper(), 'N')
        reverse_complement.append(comp_base) 
        
    reverse_complement.reverse()
    reverse_complement = ''.join(reverse_complement)
    return reverse_complement 

# function to convert the fasta file into a string containing DNA sequence 
def load_fasta(file_path):

    record = next(SeqIO.parse(file_path, "fasta")) # parse sends the file to the RAM chunk by chunk using an iterator
    # we use next because a single fasta file contains multiple genomes, like a main chromosome and maybe some "plasmids"
    # next gives us the very first block of DNA it finds 
    # record now holds a complex Biopython object that contains header name , description, DNA sequence, etc

    clean_sequence = str(record.seq).upper() 
    # we use seq to get only the DNA , the convert to upper case string , as record.seq returned a Seq object

    return clean_sequence

sequence = load_fasta("ecoli_reference.fasta") # load the fasta file 

# function to build the map from k-mer to its locations in the reference genome
def build_k_mer_index(sequence, k):

    k_mer_index = defaultdict(list) # dictionary where each key is a k-mer , value is the positionof the k-mer in the reference genome
    for i in range(len(sequence) - k+1):

        k_mer = sequence[i: i+k]

        # skip 'N' containing k-mers as it doesnt have a 2 bit binary representation and it wont be present in the LUT or the reads
        if 'N' in k_mer:
            continue

        reverse_comp_k_mer = get_reverse_complement(k_mer)

        if k_mer <= reverse_comp_k_mer: # pick between k-mer and its reverse complement
            k_mer_index[k_mer].append(i)
        else:
            k_mer_index[reverse_comp_k_mer].append(i)

    return k_mer_index

k = 11

k_mer_index = build_k_mer_index(sequence, k)

# calculate some stats to verify and also store the k-mer index in a file

# calc stats
total_unique = len(k_mer_index) # each key each unique k-mer so length gives no of unique k-mers
frequencies = [len(pos_list) for pos_list in k_mer_index.values()] # we count the number of locations in each k-mers value using pos_list,
# which helps us iterate through the values of the k-mer index . k-mer index.values gives values of the dictionary , we need the length of each one 

max_freq = max(frequencies)
mean_freq = sum(frequencies) / total_unique 

# saving to disk 
''' The pickle library performs Serialization. It takes your complex, 3D Python Dictionary living in RAM and 
crushes it down into a raw stream of 1s and 0s that can safely be written to a hard drive. Later, when you load the .pkl file, 
Python will "un-pickle" those bytes and instantly reconstruct the exact dictionary structure back into RAM! '''

output_filename = f"k_mer_index_k_{11}.pkl"

standard_dict = dict(k_mer_index) # we intially built using defaultdict, which builds data fast but contains extra metadata,
# we convert it to a standard dict to save space

with open(output_filename, 'wb') as f: # opens file stream to hard drive 
    # wb flag means write binary 
    pickle.dump(standard_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    # dump forces dictionary through serialization and writes into a file stream using f( file object or file handle)
    # highest protocol tells python to use most modern compression algo available, so we take min space 

print(f"Total unique k-mers: {total_unique}")
print(f"Max frequency: {max_freq}")
print(f"Mean frequency: {mean_freq}")








