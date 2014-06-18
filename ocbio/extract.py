# Feature extraction code

import os, time, subprocess, csv, shelve, re, pdb

class FeatureVectorAssembler():
    '''Assembles feature vectors from protein pair files, data source lists and gold standard protein pair lists.'''
    def __init__(self,sourcetab):
        #instantiate protein pair parsers
        # first parse the data source table
        # store the directory of the table and it's name
        self.sourcetabdir,self.tabfile = os.path.split(sourcetab)
        # open the table and parse for initialisation options
        c = csv.reader(open(sourcetab), delimiter="\t")
        # iterate over lines adding to list of protein pair parsers
        self.parserinitlist = []
        for line in c:
            #store the information in a dictionary
            d = {}
            d["data path"] = os.path.join(self.sourcetabdir,line[0])
            d["output path"] = os.path.join(self.sourcetabdir,line[1])
            #store options in a dictionary in the dictionary
            # if there are options
            d["options"] = {}
            if line[2] != "":
                options = line[2].split(";")
                for x in options:
                    #split each option to find out which option it is:
                    x = x.split("=")
                    #store it in the dictionary
                    # if there are invalid options this code WILL NOT DETECT THEM
                    d["options"][x[0]]= x[1]
                #update the script directory
                if "script" in d["options"].keys():
                    d["options"]["script"] = os.path.join(self.sourcetabdir,d["options"]["script"])
                #parse protindexes and validexes:
                if "protindexes" in d["options"].keys():
                    d["options"]["protindexes"] = tuple(int(v) for v in re.findall("[0-9]+", d["options"]["protindexes"]))
                if "valindexes" in d["options"].keys():
                    d["options"]["valindexes"] = tuple(int(v) for v in re.findall("[0-9]+", d["options"]["valindexes"]))
            #copy the dictionary into the list
            self.parserinitlist.append(d.copy())
        #then initialise each of these parsers and keep them in a list
        self.parserlist = []
        for parser in self.parserinitlist:
            self.parserlist.append(ProteinPairParser(parser["data path"],
                                                     parser["output path"],
                                                     **parser["options"]))
        return None

    def regenerate(self, force=False):
        '''Calls all known protein parsers and gets them to regenerate their output, if they have to.'''
        for parser in self.parserlist:
            parser.regenerate(force)
        return None
    
    def assemble(self, pairfile, outputfile, pairlabels=False, delim="\t", missinglabel="missing"):
        '''Assembles a file of feature vectors for each protein pair in a protein pair file supplied.
        
        Assumes the pairfile is tab delimited.'''
        # first parse the pairfile into a list of frozensets
        pairs = map(lambda l: frozenset(l),csv.reader(open(pairfile), delimiter="\t"))
        # open the file to put the feature vector in
        c = csv.writer(open(outputfile, "w"), delimiter=delim)
        #open all the databases and put them in a dictionary
        dbdict = {}
        for parser in self.parserinitlist:
            dbdict[parser["output path"]] = openpairshelf(parser["output path"])
        
        # then iterate through the pairs, querying all parser databases and building a list of rows
        rows = []
        for pair in pairs:
            row = []
            if pairlabels==True:
                lpair = list(pair)
                if len(lpair)==1:
                    lpair = lpair*2
                row = row+lpair
            for parser in self.parserinitlist:
                #if there are features there then append them to the row
                try:
                    row.append(dbdict[parser["output path"]][pair])
                except KeyError:
                    row.append([missinglabel])
            #copy the row into the list of rows
            rows.append(row[:])
        #check if lengths of missing labels is correct
        columns = zip(*rows)
        for col,colindex in zip(columns,range(len(columns))):

        #write all the rows to file
        c.writerows(rows)
            
        #close all the databases
        for parser in self.parserinitlist:
            dbdict[parser["output path"]].close()
        
        return None


class ProteinPairDB(shelve.DbfilenameShelf):
    '''A simple database for protein pairs using shelve.'''
    def __setitem__(self,key,value):
        #key will be frozenset so make it a list first
        key = list(key)
        #then make it a string
        if len(key) == 1:
            key = key[0]*2
        else:
            key = key[0] + "\t" + key[1]
        shelve.DbfilenameShelf.__setitem__(self,key,value)
        return None
    
    def __getitem__(self,key):
        #make two strings from the key
        key = list(key)
        if len(key) == 1:
            key1 = key[0]*2
            key1 = key[0]*2
        else:
            key1 = key[0] + "\t" + key[1]
            key2 = key[1] + "\t" + key[0]
        #try the first one
        try:
            value = shelve.DbfilenameShelf.__getitem__(self,key1)
        except KeyError:
            #couldn't find the first key, try the second
            value = shelve.DbfilenameShelf.__getitem__(self,key2)
            #if we don't find this one then error out as usual
        return value


class ProteinPairParser():
    '''Does simple parsing on data files to produce protein pair files with feature values'''
    def __init__(self,datadir,outdir,protindexes=(0,1),valindexes=[2],script=None,csvdelim="\t",ignoreheader=0):
        #first, store the initialisation
        self.datadir = datadir
        self.outdir = outdir
        self.protindexes=protindexes
        # had to hack this together from the list above
        # passing tuple in as default did not work
        self.valindexes=tuple(valindexes)
        self.script=script
        self.csvdelim=csvdelim
        self.ignoreheader=ignoreheader
        return None
    
    def regenerate(self, force=False):
        '''Regenerate the pair file from the data source
        if the data source is newer than the pair file'''
        # so first check the ages of both files
        datamtime = os.stat(self.datadir)[-2]
        if os.path.isfile(self.outdir):
            pairmtime = os.stat(self.outdir)[-2]
        else:
            #bit of a hack
            pairmtime = 0
        #if the data modification time is greater than output modification time
        if datamtime > pairmtime or force==True:
            # now regenerate the data file according to the options defined above:
            print "data file is newer than pair file"
            #if there's a script to run
            if self.script != None:
                #then execute the script
                retcode=subprocess.call("python2 %s"%self.script, shell=True)
            #first delete out of date file, if it's there
            if os.path.isfile(self.outdir):
                os.remove(self.outdir)
            #perform simple parsing to make a file of just protein pairs and the value we care about
            #and save these using shelve
            db = openpairshelf(self.outdir)
            #open the data file
            c = csv.reader(open(self.datadir), delimiter=self.csvdelim)
            #if the header should be ignored then ignore it
            if self.ignoreheader=="1":
                print "ignoring header"
                c.next()
            for line in c:
                #each line use the protein pair as a key
                #by formatting it as a frozenset
                pair = frozenset([line[self.protindexes[0]],line[self.protindexes[1]]])
                #and the value is indexed by valindexes
                db[pair] = []
                try:
                    for i in self.valindexes:
                        db[pair].append(line[i])
                except TypeError:
                    pdb.set_trace()
            db.close()
        return None

def openpairshelf(filename, flag='c', protocol=None, writeback=False):
    """Returns a ProteinPairDB object, with similar functionality to shelve.open()"""
    return ProteinPairDB(filename, flag, protocol, writeback)
