import numpy as np

class GridCell(object):

    def __init__(self,id_):
        self.id = id_;
        self.index = np.zeros((2),dtype = np.int32);
        self.lat = None;
        self.lon = None;
        ##average pair-wise correlation between grid cell
                         ## and grid cells in domain
        self.score = -1;
        ##numpy array with dimensions equal to the dimensions of the map of your data
        ##an entry is equal to 1 if a grid cell belongs to the domain
        self.map = None;


    def __eq__(self,other):
        if(isinstance(other,Domain)):
            return self.id == other.id;
        return NotImplemented;

    def __ne__(self,other):
        x = self.__eq__(other);
        if(x is not NotImplemented):
            return not x;
        return NotImplemented

    def __hash__(self):
        return hash(self.id);

