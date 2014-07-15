####### Attribution ###############################
# This code is adapted from code used in the ENTS #
# paper, found at:                                #
#     http://ents.as.wvu.edu/index.php            #
# This code was presumably written by the         #
# authors of the paper:                           #
#     * Eli Rodgers-Melnick                       #
#     * Mark Culp                                 #
#     * Stephen P. DiFazio                        #
###################################################

import pdb

class ENTSfeatures():
    def __init__(self, subcell_info, odds_dict, proteins, domine_interactions,
                 protein_domain_dict, protein_gene_dict, entreztoensembl,
                 domain_cols, svm_subcell_cols, svm_detail_cols):
        """Object to return a feature vector for a given protein pair in Entrez format"""
        #store all required data
        #for information about how this is generated 
        #see the "Inspecting ENTS code" notebook in the opencast-bio repository
        self.subcell_info = subcell_info
        self.odds_dict = odds_dict
        self.proteins = proteins
        self.domine_interactions = domine_interactions
        self.protein_domain_dict = protein_domain_dict
        self.protein_gene_dict = protein_gene_dict
        self.entreztoensembl = entreztoensembl
        self.domain_cols = domain_cols
        self.svm_subcell_cols = svm_subcell_cols
        self.svm_detail_cols = svm_detail_cols
        return None
    
    def getfeaturevector(self, protein1, protein2):
        """prototype function to retreive ENTS feature vectors for a given protein pair"""
        keys = self.subcell_info.keys()
        # Configure svm columns for organism
        gene = keys[0]
        delete_svm_pred_cols = [x for x in self.svm_subcell_cols if x not in self.subcell_info[gene]['predictions'].keys()]
        delete_svm_detail_cols = [x for x in self.svm_detail_cols if x not in self.subcell_info[gene]['svm_info'].keys()]
        for col in delete_svm_pred_cols: svm_subcell_cols.remove(col)
        for col in delete_svm_detail_cols: svm_detail_cols.remove(col)
        # Make the true interactions part of the file
        try:
            domains = self.makeDomainString(self.protein_domain_dict[protein1], self.protein_domain_dict[protein2])
        except:
            key1 = [x for x in self.protein_domain_dict.keys() if protein1.startswith(x.split('.')[0])]
            key2 = [x for x in self.protein_domain_dict.keys() if protein2.startswith(x.split('.')[0])]
            if len(key1) == 0 or len(key2) == 0:
                return None
            else:
                domains = self.makeDomainString(self.protein_domain_dict[key1[0]], self.protein_domain_dict[key2[0]])
        # Get the subcellular localization part of the string
        try: subcells = self.makeSubcellularDict(protein1, protein2)
        except KeyError: 
            key1 = [x for x in self.protein_domain_dict.keys() if protein1.startswith(x.split('.')[0])]
            key2 = [x for x in self.protein_domain_dict.keys() if protein2.startswith(x.split('.')[0])]
            if len(key1) == 0 or len(key2) == 0:
                return None
            else:
                try:
                    subcells = self.makeSubcellularDict(key1[0],key2[0])
                except:
                    #just in case
                    return None
        return domains + subcells

    def makeSubcellularDict(self,protein1, protein2):
        gene1 = self.protein_gene_dict[protein1]
        gene2 = self.protein_gene_dict[protein2]
        svm_line1 = [self.subcell_info[gene1]['predictions'][k] for k in self.svm_subcell_cols]
        svm_line2 = [self.subcell_info[gene2]['predictions'][k] for k in self.svm_subcell_cols]
        svm_line1 += [self.subcell_info[gene1]['svm_info'][k] for k in self.svm_detail_cols]
        svm_line2 += [self.subcell_info[gene2]['svm_info'][k] for k in self.svm_detail_cols]
        return svm_line1 + svm_line2

    def getTwoListCombos(self, list1, list2):
	returnList = []
	for i in xrange(len(list1)):
		for j in xrange(len(list2)):
			if sorted([list1[i],list2[j]]) in returnList: continue
			else: returnList.append(sorted([list1[i],list2[j]]))
	return returnList

    def makeDomainString(self, domains1, domains2):
        domain_dict = {}
        # Get all potential interactions
        if len(domains1) > 0 and len(domains2) > 0:
            potential_domain_pairs = self.getTwoListCombos(domains1,domains2)
        else: potential_domain_pairs = []
        # Get the domine information
        domine_pairs = list(set([tuple(sorted(x)) for x in potential_domain_pairs if tuple(sorted(x)) in self.domine_interactions]))
        domain_dict['n_domine_pairs'] = len(domine_pairs)
        domain_dict['highest_domine_conf'] = '0'
        for pair in domine_pairs:
            if self.domine_interactions[pair] == 'HC': domain_dict['highest_domine_conf'] = 'HC'
            elif self.domine_interactions[pair] == 'MC' and domain_dict['highest_domine_conf'] != 'HC':
                domain_dict['highest_domine_conf'] = 'MC'
            elif self.domine_interactions[pair] == 'LC' and domain_dict['highest_domine_conf'] not in ['HC','MC']:
                domain_dict['highest_domine_conf'] = 'LC'
        # Get the odds information
        ############## TEST ##################
        domain_dict['lowest_odds'] = 0.
        domain_dict['not_observed'] = 0
        domain_dict['not_observed_frac'] = 1.
        ################################
        domain_dict['sum_odds'] = 0.
        domain_dict['highest_odds'] = 0.
        domain_dict['n_odds_pairs'] = 0
        for pair in potential_domain_pairs:
            pair = tuple(list(sorted(pair)))
            # Check if in the odds dictionary. If so, update domain variables
            if pair in self.odds_dict:
                domain_dict['sum_odds'] += self.odds_dict[pair]
                if self.odds_dict[pair] > domain_dict['highest_odds']:
                    domain_dict['highest_odds'] = self.odds_dict[pair]
                domain_dict['n_odds_pairs'] += 1
        ################### TEST ################
                if self.odds_dict[pair] < domain_dict['lowest_odds']: domain_dict['lowest_odds'] = self.odds_dict[pair]
            else:
                domain_dict['not_observed'] += 1
        if domain_dict['n_odds_pairs'] + domain_dict['not_observed'] > 0: 
            domain_dict['not_observed_frac'] = float(domain_dict['not_observed']) / (domain_dict['n_odds_pairs'] + domain_dict['not_observed'])
            #########################################
        return [domain_dict[k] for k in self.domain_cols]

    def highestdomineconf1ofk(self,fvector):
        """Replaces string encoding of highest_domine_conf
        with 1-of-k encoding"""
        if fvector[7] == '0':
            fvector = fvector[:7] + [1,0,0,0] + fvector[8:]
        elif fvector[7] == 'LC':
            fvector = fvector[:7] + [0,1,0,0] + fvector[8:]
        elif fvector[7] == 'MC':
            fvector = fvector[:7] + [0,0,1,0] + fvector[8:]
        elif fvector[7] == 'HC':
            fvector = fvector[:7] + [0,0,0,1] + fvector[8:]
        return fvector

    def __getitem__(self,key):
        """Where key is a protein pair, returns a feature vector"""
        #the pair should be a frozenset of Entrez IDs
        pair = list(key)
        if len(pair) == 1:
            pair = pair*2
        #convert both to ensembl
        for p in self.entreztoensembl[pair[0]]:
            if p in self.proteins:
                protein1 = p
                break
            else:
                raise KeyError("Unknown protein")
        for p in self.entreztoensembl[pair[1]]:
            if p in self.proteins:
                protein2 = p
                break
            else:
                raise KeyError("Unknown protein")
        #retrieve the feature vector
        fvector = self.getfeaturevector(protein1,protein2)
        if fvector == None:
            raise KeyError("No feature vector found for pair {0}".format(pair))
        #transform feature vector
        # fix string encoding of highest_domine_conf
        fvector = self.highestdomineconf1ofk(fvector)

        return fvector
