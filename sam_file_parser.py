# this script reads a SAM file and extracts necessary info and add it into a TSV file
# TSV stands for Tab Separated Values, its a simple text file where each line is a row and columns are separated by tabs

def parse_sam_truth(sam_file, output_tsv):

    # open both the input SAM file and the output TSV file at the same time
    with open(sam_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        
        # write the header row for our new TSV file so we know what the columns are
        outfile.write("read_id\tref_start\tref_end\tstrand\tread_length\n")

        # read the SAM file line-by-line, to avoid loading entire file into memory at once, which will be a problem for large files
        for line in infile:
            
            # skip header lines(they all start with '@')
            if line.startswith('@'):
                continue
                
            # clean up hidden newline characters and split by the Tab character
            columns = line.strip().split('\t') # strip removes hidden newline characters, and then we split by tab
            
            # make sure the line actually has enough data
            if len(columns) < 10:
                continue

            # extract columns
            read_id = columns[0]
            flag = int(columns[1]) # convert to int for math later
            pos = int(columns[3])  
            seq = columns[9]

            # do some calculations to get start and end position of read in reference genome
            ref_start = pos - 1 # sam file stores positions in 1-based indexing, but we want 0-based indexing for further processing
            read_length = len(seq)
            ref_end = ref_start + read_length

            # flag value is 16 if the read came from the reverse strand
            if flag & 16:
                strand = '-' # reverse strand
            else:
                strand = '+' # forward strand

            # write the calculated data to our new TSV file, separated by tabs
            outfile.write(f"{read_id}\t{ref_start}\t{ref_end}\t{strand}\t{read_length}\n")


parse_sam_truth("02_simulated_short.sam", "02_short_read_truth.tsv")
print("Finished parsing SAM file!")