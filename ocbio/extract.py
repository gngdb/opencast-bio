# Feature extraction code

import os
import time
import subprocess
import csv
import shelve
import re
import sys


class FeatureVectorAssembler():
    '''Assembles feature vectors from protein pair files, data source lists
    and gold standard protein pair lists.'''
    def __init__(self, sourcetab, verbose=False):
        # Instantiate protein pair parsers
        # first parse the data source table
        # store the directory of the table and it's name
        self.sourcetabdir, self.tabfile = os.path.split(sourcetab)

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
                # parse protindexes and validexes:
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
                                                     **parser["options"]))
        v_print("Finished Initialisation.")
        return None

    def regenerate(self, force=False, verbose=False):
        '''Calls all known protein parsers and gets them to regenerate their
        output, if they have to.'''
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
        v_print("Reading pairfile: {0}".format(pairfile))
        # first parse the pairfile into a list of frozensets
        pairs = map(lambda l: frozenset(l), csv.reader(open(pairfile), delimiter="\t"))
        # open the file to put the feature vector in
        c = csv.writer(open(outputfile, "w"), delimiter=delim)
        # open all the databases and put them in a dictionary
        dbdict = {}

        v_print("Opening databases:")
        for parser in self.parserinitlist:
            dbdict[parser["output path"]] = openpairshelf(parser["output path"])
            v_print("\t {0} open".format(parser["output path"]))

        v_print("Checking feature sizes:")
        # check size of feature in each file
        # will be important later
        featuresizes = {}
        for parser in self.parserinitlist:
            k = dbdict[parser["output path"]].keys()[0]
            featuresizes[parser["output path"]] = len(dbdict[parser["output path"]][k])
            v_print("\t Database {0} contains features of size {1}.".format(parser["output path"],
                                                                            featuresizes[parser["output path"]]))
        if verbose:
            sys.stdout.write("Writing feature vectors")
            lcount = 0
            # counters for each database reporting numbers of missing values
            mcount = {}
            for parser in self.parserinitlist:
                mcount[parser["output path"]] = 0
        # then iterate through the pairs, querying all parser databases and building a list of rows
        rows = []
        for pair in pairs:
            row = []
            if pairlabels is True:
                lpair = list(pair)
                if len(lpair) == 1:
                    lpair = lpair * 2
                row = row + lpair
            for parser in self.parserinitlist:
                # if there are features there then append them to the row
                try:
                    row = row + dbdict[parser["output path"]][pair]
                except KeyError:
                    row = row + [missinglabel] * featuresizes[parser["output path"]]
                    if verbose:
                        mcount[parser["output path"]] += 1
            c.writerow(row)
            if verbose:
                lcount = lcount+1
                if lcount % 1000 == 0:
                    sys.stdout.write(".")
        if verbose:
            sys.stdout.write("\n")
            print "Wrote {0} vectors.".format(lcount)
            for parser in self.parserinitlist:
                percentage_match = 100.0 - 100.0 * mcount[parser["output path"]] / lcount
                print "Matched {0:.2f} % of protein pairs in {1} to {2}".format(percentage_match,
                                                                                pairfile,
                                                                                parser["output path"])
        # close all the databases
        for parser in self.parserinitlist:
            dbdict[parser["output path"]].close()

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
            key1 = key[0] * 2
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
                 ignoreheader=0):
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
        return None

    def regenerate(self, force=False, verbose=False):
        '''Regenerate the pair file from the data source
        if the data source is newer than the pair file'''
        # so first check the ages of both files
        datamtime = os.stat(self.datadir)[-2]
        if os.path.isfile(self.outdir):
            pairmtime = os.stat(self.outdir)[-2]
        else:
            # bit of a hack
            pairmtime = 0
        # if the data modification time is greater than output modification time
        if datamtime > pairmtime or force is True:
            # now regenerate the data file according to the options defined above:
            if verbose and datamtime > pairmtime:
                if pairmtime == 0:
                    print "Database file not found, regenerating at {0} from {1}.".format(self.outdir, self.datadir)
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
            # first delete out of date file, if it's there
            if os.path.isfile(self.outdir):
                os.remove(self.outdir)
            # perform simple parsing to make a file of just protein pairs and the value we care about
            # and save these using shelve
            db = openpairshelf(self.outdir)
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

                db[pair] = values[:]

                if verbose:
                    lcount = lcount + 1
                    if lcount % 1000 == 0:
                        sys.stdout.write(".")

            if verbose:
                sys.stdout.write("\n")
                print "Parsed {0} lines.".format(lcount)
            db.close()

        return None


def openpairshelf(filename, flag='c', protocol=None, writeback=False):
    """Returns a ProteinPairDB object, with similar functionality to shelve.open()"""
    return ProteinPairDB(filename, flag, protocol, writeback)
