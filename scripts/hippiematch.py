#!/usr/bin/python2

import sys
import csv
import os

def main(hippiefile, pairfile, out):
    """A script to take a file output by the HIPPIE script and match the pairs in it to a file containing protein pairs:
    input:
        hippiefile - output of HIPPIE script
        pairfile - list of protein pairs
    output:
        out - all the pairs it could find in the HIPPIE file with their confidence values"""
    # Reading in the HIPPIE file
    #initialise csv reader
    c = csv.reader(open(hippiefile), delimiter="\t")
    #make dictionary using frozensets as keys with the confidence scores as values
    hippieids = {}
    for line in c:
        k = frozenset([line[1],line[3]])
        hippieids[k] = line[4]

    # Reading in the file of protein pairs
    #initialise csv reader
    c = csv.reader(open(pairfile), delimiter="\t")
    #make dictionary using frozensets as keys:
    posids = {}
    for line in c:
        line = frozenset(line)
        posids[line] = 1

    # Write new file:
    #then rewrite the training.positive.HIPPIE.txt file:
    c = csv.writer(open(out, "w"), delimiter="\t")
    for k in hippieids.keys():
        try:
            if posids[k]:
                l = list(k)
                try:
                    c.writerow([l[0],l[1],hippieids[k]])
                except:
                    #ignore self-interactions
                    pass
        except KeyError:
            #ignore missing pairs
            pass

    #how many lines does the new file have?
    count = int(os.popen("wc -l < " + out).read())
    # Report how well it went:
    print "%i of %i pairs matched."%(count,len(posids.keys()))

    return None


if __name__=="__main__":
    if sys.argv[1] == "-h":
        print "Usage: python2 hippiematch.py hippiefile pairfile outputfile"
        print "Where: "
        print "    hippiefile is the out of the HIPPIE script"
        print "    pairfile is the file of protein pairs to match"
        print "    outputfile is the name of the file to output to"
    else:
        main(*sys.argv[1:])
