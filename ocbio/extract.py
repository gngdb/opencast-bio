# Feature extraction code
# Header pending

import os
import time
import subprocess
import csv
import shelve
import re
import sys
import pickle
import geneontology
#import pdb
import time
import ents
import ppipred
import irefindex

def verbosecheck(verbose):
    '''returns a function depending on the state of the verbose flag'''
    if verbose:
        def v_print(*args):
            '''declare v_print function that prints to stdout
            if verbose flag is on'''
            for argument in args:
                print argument,
                print
    else:
        def v_print(*args):
            None
    return v_print

class FeatureVectorAssembler():
    '''Assembles feature vectors from protein pair files, data source lists
    and gold standard protein pair lists.'''
    def __init__(self, sourcetab, verbose=False):
        # Instantiate protein pair parsers
        # first parse the data source table
        # store the directory of the table and it's name
        self.sourcetabdir, self.tabfile = os.path.split(sourcetab)

        v_print = verbosecheck(verbose)

        v_print("Using {0} from top data directory {1}.".format(self.sourcetabdir,
                                                                self.tabfile))

        # open the table and parse for initialisation options
        c = csv.reader(open(sourcetab), delimiter="\t")
        # iterate over lines adding to list of protein pair parsers
        v_print("Reading data source table:")

        self.parserinitlist = []
        for line in c:
            # store the information in a dictionary
            d = {}
            d["data path"] = os.path.join(self.sourcetabdir, line[0])
            d["output path"] = os.path.join(self.sourcetabdir, line[1])
            # store options in a dictionary in the dictionary
            # if there are options
            d["options"] = {}
            if line[2] != "":
                options = line[2].split(";")
                for x in options:
                    # split each option to find out which option it is:
                    x = x.split("=")
                    # store it in the dictionary
                    # if there are invalid options this code
                    # WILL NOT DETECT THEM
                    d["options"][x[0]] = x[1]
                # update the script directory
                if "script" in d["options"].keys():
                    d["options"]["script"] = os.path.join(self.sourcetabdir,
                                                          d["options"]["script"])
                # parse protindexes and valindexes:
                if "protindexes" in d["options"].keys():
                    d["options"]["protindexes"] = tuple(int(v) for v in re.findall("[0-9]+", d["options"]["protindexes"]))
                if "valindexes" in d["options"].keys():
                    d["options"]["valindexes"] = tuple(int(v) for v in re.findall("[0-9]+", d["options"]["valindexes"]))
            # copy the dictionary into the list

            v_print("\t"+"Data source: {0} to be processed to {1}".format(d["data path"],
                                                                          d["output path"]))

            self.parserinitlist.append(d.copy())
        # then initialise each of these parsers and keep them in a list
        self.parserlist = []
        v_print("Initialising parsers.")
        for parser in self.parserinitlist:
            self.parserlist.append(ProteinPairParser(parser["data path"],
                                                     parser["output path"],
                                                     verbose=verbose,
                                                     **parser["options"]))
        v_print("Finished Initialisation.")
        return None

    def regenerate(self, force=False, verbose=False):
        '''Calls all known protein parsers and gets them to regenerate their
        output, if they have to.'''
        v_print = verbosecheck(verbose)
        v_print("Regenerating parsers:")
        for parser in self.parserlist:
            v_print("\t parser {0}".format(self.parserlist.index(parser)))
            parser.regenerate(force, verbose)
        return None

    def assemble(self, pairfile, outputfile, pairlabels=False, delim="\t",
                 missinglabel="missing", verbose=False):
        '''Assembles a file of feature vectors for each protein pair in a
        protein pair file supplied.
        Assumes the pairfile is tab delimited.'''
        v_print = verbosecheck(verbose)
        v_print("Reading pairfile: {0}".format(pairfile))
        # first parse the pairfile into a list of frozensets
        pairs = map(lambda l: frozenset(l), csv.reader(open(pairfile), delimiter="\t"))
        # open the file to put the feature vector in
        c = csv.writer(open(outputfile, "w"), delimiter=delim)

        v_print("Checking feature sizes:")
        # check size of feature in each file
        # will be important later
        featuresizes = {}
        for parser, i in zip(self.parserlist, range(len(self.parserlist))):
            #try to get an example feature
            examplefeature = None
            for pair in pairs:
                try:
                    examplefeature = parser[pair]
                except KeyError:
                    #keep trying
                    pass
                #if examplefeature != None:
                #    break
            #check we actually got an example feature
            if examplefeature == None:
                # should probably not include a feature that's going to be all missing values
                del self.parserlist[i]
                v_print("\t Feature from {0} does not map to these protein pairs.".format(parser.datadir))
            else:
                #then we've got a feature so we should see what size it is
                featuresizes[parser.datadir] = len(examplefeature)
                v_print("\t Data source {0} produces features of size {1}.".format(parser.datadir,
                                                                            featuresizes[parser.datadir]))
        if verbose:
            sys.stdout.write("Writing feature vectors")
            lcount = 0
            # counters for each database reporting numbers of missing values
            mcount = {}
            for parser in self.parserlist:
                mcount[parser.datadir] = 0
        # then iterate through the pairs, querying all parser databases and building a list of rows
        rows = []
        for pair in pairs:
            row = []
            if pairlabels is True:
                lpair = list(pair)
                if len(lpair) == 1:
                    lpair = lpair * 2
                row = row + lpair
            for parser in self.parserlist:
                # if there are features there then append them to the row
                try:
                    row = row + parser[pair]
                except KeyError:
                    row = row + [missinglabel] * featuresizes[parser.datadir]
                    if verbose:
                        mcount[parser.datadir] += 1
            c.writerow(row)
            if verbose:
                lcount = lcount+1
                if lcount % 10000 == 0:
                    sys.stdout.write(".")
        if verbose:
            sys.stdout.write("\n")
            print "Wrote {0} vectors.".format(lcount)
            for parser in self.parserlist:
                percentage_match = 100.0 - 100.0 * mcount[parser.datadir] / lcount
                print "Matched {0:.2f} % of protein pairs in {1} to features from {2}".format(percentage_match,
                                                                                pairfile,
                                                                                parser.datadir)
        return None

    def close(self, verbose=False):
        v_print = verbosecheck(verbose)
        for parser in self.parserlist:
            if parser.db != None:
                parser.close()
                v_print("{0} closed".format(parser.outdir))
        return None


class ProteinPairDB(shelve.DbfilenameShelf):
    '''A simple database for protein pairs using shelve.'''
    def __setitem__(self, key, value):
        # key will be frozenset so make it a list first
        key = list(key)
        # then make it a string
        if len(key) == 1:
            key = key[0] * 2
        else:
            key = key[0] + "\t" + key[1]
        shelve.DbfilenameShelf.__setitem__(self, key, value)
        return None

    def __getitem__(self, key):
        # make two strings from the key
        key = list(key)
        if len(key) == 1:
            key1 = key[0] * 2
            key2 = key[0] * 2
        else:
            key1 = key[0] + "\t" + key[1]
            key2 = key[1] + "\t" + key[0]
        # try the first one
        try:
            value = shelve.DbfilenameShelf.__getitem__(self, key1)
        except KeyError:
            # couldn't find the first key, try the second
            value = shelve.DbfilenameShelf.__getitem__(self, key2)
            # if we don't find this one then error out as usual
        return value

    def keys(self):
        # retrieve the string keys used by shelve
        ks = shelve.DbfilenameShelf.keys(self)
        # convert them to frozensets
        ks = map(lambda x: frozenset(x.split("\t")), ks)
        return ks


class ProteinPairParser():
    '''Does simple parsing on data files to produce protein pair files with feature values'''
    def __init__(self,
                 datadir,
                 outdir,
                 protindexes=(0, 1),
                 valindexes=[2],
                 script=None,
                 csvdelim="\t",
                 ignoreheader=0,
                 generator=False,
                 verbose=False):
        v_print = verbosecheck(verbose)
        # first, store the initialisation
        self.datadir = datadir
        self.outdir = outdir
        self.protindexes = protindexes
        # had to hack this together from the list above
        # passing tuple in as default did not work
        self.valindexes = tuple(valindexes)
        self.script = script
        self.csvdelim = csvdelim
        self.ignoreheader = ignoreheader
        if generator:
            #then open up this pickle file
            f = open(generator)
            self.generator = pickle.load(f)
            f.close()
            self.db = None
        else:
            #otherwise open the database that is assumed to exist 
            self.generator = None
            self.db = openpairshelf(self.outdir)
            #check if this is a new database
            try:
                v_print("Database {0} last updated {1}".format(self.outdir,time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(self.db["last written"]))))
            except KeyError:
                v_print("Database {0} may be empty, must be updated, please run regenerate.".format(self.outdir))
        return None

    def regenerate(self, force=False, verbose=False):
        '''Regenerate the pair file from the data source
        if the data source is newer than the pair file'''
        v_print = verbosecheck(verbose)
        if self.generator == None:
            # so first check the age of the data file
            datamtime = os.stat(self.datadir)[-2]
            # check if the database file has ever been written to before
            try:
                v_print("Database {0} last updated {1}".format(self.outdir,time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(self.db["last written"]))))
            except KeyError:
                #if not make sure it goes through the next if statement
                self.db["last written"] = 0
            # if the data modification time is greater than last time we wrote to the database
            if datamtime > self.db["last written"] or force is True:
                # now regenerate the data file according to the options defined above:
                if verbose and datamtime > self.db["last written"]:
                    if self.db["last written"] == 0:
                        print "Database file may be empty, regenerating at {0} from {1}.".format(self.outdir, self.datadir)
                    else:
                        print "Data file {0} is newer than processed database {1}, regenerating.".format(self.datadir, self.outdir)
                if verbose and force:
                    print "Forcing regeneration of database {0} from data file {1}.".format(self.outdir, self.datadir)
                # if there's a script to run
                if self.script is not None:
                    v_print("Executing script: {0}.".format(self.script))
                    # then execute the script
                    retcode = subprocess.call("python2 {0}".format(self.script), shell=True)

                    v_print("Script returned: {0}".format(retcode))
                # open the data file
                c = csv.reader(open(self.datadir), delimiter=self.csvdelim)
                # if the header should be ignored then ignore it

                if self.ignoreheader == "1":
                    v_print("Ignoring header.")
                    c.next()

                if verbose:
                    sys.stdout.write("Filling database")
                    lcount = 0

                for line in c:
                    # each line use the protein pair as a key
                    # by formatting it as a frozenset
                    pair = frozenset([line[self.protindexes[0]], line[self.protindexes[1]]])
                    # and the value is indexed by valindexes
                    values = []

                    for i in self.valindexes:
                        values.append(line[i])

                    self.db[pair] = values[:]

                    if verbose:
                        lcount = lcount + 1
                        if lcount % 1000 == 0:
                            sys.stdout.write(".")

                if verbose:
                    sys.stdout.write("\n")
                    print "Parsed {0} lines.".format(lcount)
                    # add the current time to the database "last written" entry:
                    self.db["last written"] = time.time()
        else:
            v_print("Custom generator function, no database to regenerate.")
        return None

    def __getitem__(self,key):
        if self.generator != None:
            #try and read a key from the custom generator
            return self.generator[key]
        else:
            #read key from database
            return self.db[key]

    def close(self):
        self.db.close()
        return None


def openpairshelf(filename, flag='c', protocol=None, writeback=False):
    """Returns a ProteinPairDB object, with similar functionality to shelve.open()"""
    return ProteinPairDB(filename, flag, protocol, writeback)
