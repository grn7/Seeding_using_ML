# this file is for parsing the MAF file to extract info for long reads and write it to a TSV file 

def parse_maf_truth(maf_file, output_tsv):
    # open the input MAF and output TSV files simultaneously
    with open(maf_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        
        # write the TSV header row
        outfile.write("read_id\tref_start\tref_end\tstrand\tread_length\n")

        # set up "memory" variables(this holds our state)
        current_ref_start = None
        current_ref_length = None
        current_strand = None

        # read line by line
        for line in infile:
            
            # use split() with no arguments to automatically split by any amount of spaces/tabs (whitespaces basically)
            columns = line.strip().split()
            
            # skip empty lines
            if not columns:
                continue

            # check if it's a new alignment block
            if columns[0] == 'a':
                # RESET our memory variables for the new block!
                current_ref_start = None
                current_ref_length = None
                current_strand = None
                continue
                
            # check if it's a Sequence ('s') line 
            elif columns[0] == 's':
                
                # if memory is empty, this must be the first 's' line (the Reference Genome)
                if current_ref_start is None: 
                    # store the data in our memory variables
                    current_ref_start = int(columns[2])
                    current_ref_length = int(columns[3])
                    current_strand = columns[4]

                    # extract the full length of the E. coli chromosome (column 5)
                    src_size = int(columns[5])

                    # MAF reverse strand normalization
                    if current_strand == '-':
                        current_ref_start = src_size - current_ref_start - current_ref_length

                # if memory has data, this must be the second 's' line (the Simulated Read)
                else:
                    # extract read specific data
                    read_id = columns[1]
                    read_length = int(columns[5]) 
                    
                    # calculate where it ends on the reference
                    ref_end = current_ref_start + current_ref_length
                    
                    # write info to file
                    outfile.write(f"{read_id}\t{current_ref_start}\t{ref_end}\t{current_strand}\t{read_length}\n")


parse_maf_truth("02_simulated_long_0001.maf", "02_long_read_truth.tsv")
print("Finished parsing Long Read MAF file!")