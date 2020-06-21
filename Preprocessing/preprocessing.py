import sys
import os
import time
import numpy as np

def remove_mean(climate_field):
	"""
	Normalizes time series to zero mean
	Input:
		climate_field: Numpy mask array. The climate field time series
	Returns:
		climate_filed: Numpy mask array. The climate field zero mean time series
	"""
	mean = np.mean(climate_field,axis=0);
	climate_field = climate_field - mean;
	return climate_field;



def remove_variance(climate_field):
    """ 
    Normalizes time series to unit variance
    Input:
        climate_field: Numpy mask array. The climate field time series
    Returns:
        climate_filed: Numpy mask array. The climate field unit variance time series
    """
    std = np.std(climate_field,axis = 0);
    climate_field = climate_field/std;
    return climate_field;
