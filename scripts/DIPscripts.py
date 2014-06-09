def DIPgetproteinIDs(DIPintfile):
    """ Takes the filename of a DIP interaction file and writes three files containing the proteins IDs
    found in that DIP interaction file. Returns uniprotIDs (because you need it to run the next thing)"""
    import csv
    from collections import OrderedDict
    from matplotlib.cbook import flatten
    c = csv.reader(open(DIPintfile), delimiter="\t")
    #skip first line
    c.next()
    #take out the first two strings for every line
    IDstrings = map(lambda x: (x[0],x[1]), c )
    # The below uniqueifier is from [this stackoverflow post](http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order).
    # would quite like to flatten this list
    # then remove duplicate entries
    IDstrings = list(OrderedDict.fromkeys(flatten(IDstrings)))
    print "Number of proteins in DIP dataset is %i"%len(IDstrings)
    print "Example entries:"
    print IDstrings[0:10]
    # Next have to remove entries that only have DIP identifiers because I can't work with those.
    # Or I could just split all the strings by the `|` and split it up into the different identifier types and push all three files through the uniprot web system.
    splitstrings = map(lambda x: x.split("|"), IDstrings)
    # then have to sort into three lists
    # get ready for some super readable code
    # who needs programming efficiency?
    DIPIDs = []
    uniprotIDs = []
    refseqIDs = []
    for line in splitstrings:
        #can only look at the strings themselves to sort
        for s in line:
            #test for each case:
            if "DIP" in s:
                DIPIDs.append(s)
            elif "uniprot" in s:
                uniprotIDs.append(s)
            elif "refseq" in s:
                refseqIDs.append(s)
    #write to three files
    csv.writer(open("IDs.DIP.txt", "w"),delimiter="\n").writerow(DIPIDs)
    uniprotIDs = map(lambda x: x.split(":")[1],uniprotIDs)
    csv.writer(open("IDs.uniprot.txt", "w"),delimiter="\n").writerow(uniprotIDs)
    refseqIDs = map(lambda x: x.split(":")[1],refseqIDs)
    csv.writer(open("IDs.refseq.txt", "w"),delimiter="\n").writerow(refseqIDs)
    return uniprotIDs

def gettotaluniprot(DIPmap,refseqmap,uniprotIDs):
    """ Takes filenames for DIP and refseq map tables from uniprot and a list of uniprotIDs and combines the
    uniprotIDs in all three then writes it to a file and returns it"""
    import csv
    from collections import OrderedDict
    import pdb
    # Using the uniprot website to map these to Entrez.
    # Starting with the DIP IDs, got two files out:
    # 
    # 1. A mapping table from DIP to uniprot ID
    # 2. A list of targets
    # 
    # Interestingly, though, these two files are completely different lengths. Which could mean that a number of the DIP IDs map to the same uniprot ID. Checking if this is the case:

    # In[131]:

    #load in DIP mapping table
    c = csv.reader(open(DIPmap), delimiter="\t")
    #ignore first row
    c.next()
    # convert to lists
    DIPmap = list(c)
    # unpack and zip to get uniprot
    u = zip(*DIPmap)[1]
    # how long is it to start with
    print "Before removing duplicates length of uniprot is %i"%len(u)
    # then remove duplicates
    u = list(OrderedDict.fromkeys(u))
    # how long is it now?
    print "After removing duplicates length of uniprot is %i"%len(u)

    # Appears some number of these even have the same name so seems like it shouldn't be a syntax problem.
    # Should probably remove duplicates before sending this to uniprot.
    # Some of the entries also have extra numbers after dashes I could try slicing off as well.

    # Putting this into uniprot mapped 2,135 of 2,189.
    # Using this mapping table to add to the uniprot list as before:

    #load in refseq mapping table
    c = csv.reader(open(refseqmap), delimiter="\t")
    #ignore first row
    c.next()
    # convert to lists
    refseqmap = list(c)
    # unpack and zip to get uniprot
    ur = zip(*refseqmap)[1]
    # how long is it to start with
    print "Before removing duplicates length of uniprot is %i"%len(u)
    # then remove duplicates
    ur = list(OrderedDict.fromkeys(u))
    # how long is it now?
    print "After removing duplicates length of uniprot is %i"%len(u)


    # __Wait a second__, how is the mapping table longer than the number of identifiers mapped? Is this just the general mapping table for all proteins or something?
    # 
    # __Nope__, looks like the refseq IDs map onto multiple uniprot IDs...
    # 
    # Combining everything together into a real total uniprot list:

    # In[139]:

    #remove duplicates
    totaluniprot = list(OrderedDict.fromkeys(u+uniprotIDs+ur))
    #slice off anything extra
    totaluniprot = map(lambda x: x.split("-")[0], totaluniprot)

    print "Length before removing duplicates was %i and after is now %i."%(len(u+uniprotIDs+ur),len(totaluniprot))

    #write the file
    csv.writer(open("total.uniprot.txt", "w"),delimiter="\n").writerow(totaluniprot)

    return totaluniprot
