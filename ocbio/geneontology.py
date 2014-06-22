# coding: utf-8
import numpy as np

class GOfvectors():
    def __init__(self,goEntrezdict,terms):
        self.goEntrezdict = goEntrezdict
        self.terms = terms
        return None
    
    def __getitem__(self,key):
        pair = list(key)
        if len(pair)==1:
            pair = pair*2
        # First define the domains, I'm guessing they're not changing
        domains = ['molecular_function','cellular_component','biological_process']
        # initialise the vectors we're going to be using:
        govectordict= {}
        for entrezID in pair:
            vec = []
            for domain in domains:
                for term in self.terms[domain]:
                    try:
                        if term in self.goEntrezdict[entrezID][domain]:
                            vec.append(1)
                        else:
                            vec.append(0)
                    except KeyError:
                        vec.append(0)
            #save to dictionary
            govectordict[entrezID] = np.array(vec[:])
        #iterate over combinations of Entrez pairs
        #building the line as it goes
        line = list(govectordict[pair[0]]*govectordict[pair[1]])
        #yield this feature vector
        return line
