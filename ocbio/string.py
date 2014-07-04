# This is a class intended to store and produce feature vectors
# derived from the detailed interaction results from the
# string database as described in the notebook:
# pending link

class features():
    def __init__(self,featuredict):
        self.featuredict = featuredict
        return None

    def __getitem__(self,key):
        try:
            return self.featuredict[key]
        except KeyError:
            return ["0"]*8
        # if you can't find an entry in the dictionary return a zeroed vector
