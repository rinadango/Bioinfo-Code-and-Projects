# The file Q3-input.xml containts BLAST results for Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
# complete genome (available in the NCBI Nucleotide entry NC_045512)(change it your BLAST xml file). 
# Given the BLAST results, the script gives the geographical 
# locations from which a sample of a virus genome has been obtained that is at least 99 % identical to the query.
# Prints one location per line

from Bio import SeqIO
import Bio.Entrez as BE
import Bio.SearchIO as BSIO
import Bio.SeqFeature as BSF

# set email globally for NCBI
BE.email = 'goodmail@gmail.com'

# The whole retrieved genbank file put into a list
# We retrevied it a bit below this line of code
seq = list(SeqIO.parse("sequences.gb","gb"))
#print(seq)

# Function for percent identity that is more or equal to 99
def identity(hsp):
    return (hsp.ident_num / hsp.aln_span) * 100 >= 99

# Reading the blast results and getting the ids
with open('Q3-input.xml') as file:
    # Parse file as xml
    for result in BSIO.parse(file, 'blast-xml'):
        # Getting the BLAST ids that have 99 percent identity
        filtered_result = result.hsp_filter(identity)
        #print(filtered_result)
        
        # Empty list to store our BLAST IDs
        # Prints out all of the results that were filtered
        db_seq_ids = []
        
        # Looks trough the hits what we have filtered out and then gets only the IDS to our empty list made above this loop
        for hit in filtered_result:
            idd = hit.accession
            db_seq_ids.append(idd)
        
        # Fetch records from "nucleotide" database using our gotten IDs 'db_seq_ids'
        handle = BE.efetch(db="nucleotide",
                    id=db_seq_ids,
                    # return sequences in GenBank plain text format
                    rettype="gb",
                    retmode="text")
            
        # Save all of the results into 1 file
        with open('sequences.gb', 'w') as f:
            f.write(handle.read())


# For every BLAST ID or NCBI ID in our genbank file
for s in seq:
    print()
    
    # Name of BLAST hit IDs or NCBI IDs
    names = s.name
    print("BLAST hit IDs/NCBI IDs::",names)

    # Finds all of the countries that are in genbank file sequences with feature - 'source'
    for f in s.features:
        # Check if the feature in genbank file is 'source'
        if f.type == "source":
            # Check if qualifier is 'country', which has our country name.
            if "country" in f.qualifiers:
                # Get all of those countries
                ctr = f.qualifiers["country"]
                # Print countries line by line nicelly
                for country in ctr:
                    print("Country::",country)