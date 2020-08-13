class Domain(object):
    
    def __init__(self,domain_id):
        self.domain_id = domain_id;
        self.grid_cells = {}; ##key grid cell id, value GridCell
        self.homogeneity = -1.;
        self.map = None;
        
        
    def __eq__(self,other):
        if(isinstance(other,Domain)):
            return self.domain_id == other.domain_id;
        return NotImplemented;
    
    def __ne__(self,other):
        x = self.__eq__(other);
        if(x is not NotImplemented):
            return not x;
        return NotImplemented
    
    def __hash__(self):
        return hash(self.domain_id);
