# The function returns those amino acid residues of the protein that are post-translationally modified according to UniProt. 
# The residues are returned as a list of (position, residue) tuples where position is the zero-indexed position of the residue 
# (as integer) and residue is the single-letter code of the residue (as string).

def get_modified_residues(uid):
    # Methods for the code to work
    import requests as R
    import Bio.SeqIO as BSIO
    import collections as C
    
    # UniProt ID
    uidd = uid

    # UniProt base url
    base = "https://www.uniprot.org/uniprot/%s.%s"

    # Format
    fmt = 'xml'

    # Adds the base url to idd and xml format by replacing the %s.%s
    # It would be something like this if the id (idd) is 11:
    # https://www.uniprot.org/uniprot/11.xml
    url = base%(uidd, fmt)

    # Gets the URL response so that the contents of the URL could be saved
    response = R.get(url)

    # Prints URLs just to see how it came out
    #print(url)

    # For every id (uidd) - uniprot id we have, save the id's content into an XML file.
    # So let's say our id (idd) is 11
    # Then the xml file name would be 11.xml and the https://www.uniprot.org/uniprot/11.xml content ...
    # ... would be save into that 11.xml file
    #
    # Opens file
    output_file = open(f'{uidd}.xml', 'w')
    # Writes to file
    output_file.write(response.text)
    # Closes the file
    output_file.close()
    
    # For the xml file with the id we need
    r = BSIO.read(f'{uidd}.xml', "uniprot-xml")

    # Empty list for list of tuples
    positionRes = []

    # For every feature in our .xml file
    for feature in r.features:
        # We are looking for modified residues in the xml file
        if feature.type == "modified residue":
            #print(feature.location,feature.extract(r.seq))
            
            # feature.location gives us the location of the residues that are modified
            k = (str(feature.location),)

            # feature.extract(r.seq) gives us the residues themselves
            # We put them into a tuple: the location and the residues
            new = k + (str(feature.extract(r.seq)),)

            # Append them into a list of tuples
            positionRes.append(new)
    
    # Return the results      
    return positionRes

# input
uid = 'P60709'

# output
results = get_modified_residues(uid)
for pos, res in results:
    print("{pos:2}\t{res}".format(pos=pos, res=res))