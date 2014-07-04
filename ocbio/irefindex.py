# class for storing dictionary of feature vectors from iRefIndex
# will also ensure that zeroed vectors are returned where there is no evidence

class features()
    def __init__(self,featuredict):
        self.featuredict = featuredict
        return None

    def __getitem__(self,key):
        try:
            return self.featuredict[key]
        except KeyError:
            return ["0"]*11
        # if you can't find an entry in the dictionary return a zeroed vector
