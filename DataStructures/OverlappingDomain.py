class OverlappingDomain(object):
    
    def __init__(self, id_a, id_b):
        ##the ids of the two overlapping domains
        self.id_a = id_a;
        self.id_b = id_b;
        ##indices (x,y coordinates in map of the two overlapping domains)
        self.map = None;
        self.homogeneity = -1;
        
    def __eq__(self,other):
        if(isinstance(other,OverlappingDomain)):
            return self.id_a == other.id_a and self.id_b == other.id_b;
        return NotImplemented;
    
    def __ne__(self,other):
        x = self.__eq__(other);
        if(x is not NotImplemented):
            return not x;
        return NotImplemented
    
    def __hash__(self):
        k1 = -1;
        k2 = -1;
        if(self.id_a > self.id_b):
            k1 = self.id_a;
            k2 = self.id_b;
        else:
            k1 = self.id_b;
            k2 = self.id_a;
        pk1k2 = (k1+k2)*(k1+k2+1)/2 + k1;
        return hash(pk1k2);

