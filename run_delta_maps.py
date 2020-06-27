import sys
import os
import time
import numpy as np
import netCDF4
from netCDF4 import Dataset
from Preprocessing.preprocessing import remove_mean
from Preprocessing.preprocessing import remove_variance
from Utils import utils
from SeedIdentification.seed_identification import seed_identification
from matplotlib import pyplot as plt



def load_data(path_to_data, climate_variable, latitude_string, longitude_string):
	"""
	Loads a netcdf file containing the 3D grid
	Input:
		path_to_data: (String) path to .nc file
		climate_variable: (String) the name of the climate variable one wants to extract from the .nc file
		latitude_string: (String) the name of the latitude variable in the .nc file
		longitude_string: (String) the name of the longitude variable in the .nc file
	Returns:
		climate_field: numpy masked array of dimensions time,lat,lon
		latitudes: numpy masked array of dimension X, where X is the latitude of the dataset
		longitudes: numpy masked array of dimension X, where X is the longitude of the dataset
	"""
	nc_fid = Dataset(path_to_data, 'r');
	climate_field = nc_fid.variables[climate_variable][:];
	latitudes = nc_fid.variables[latitude_string][:];
	longitudes = nc_fid.variables[longitude_string][:];

	return climate_field, latitudes, longitudes;





if __name__ == "__main__":
	
	##convert this to user input when done
	path_to_data = "Data/COBEv2_preProcessed.nc";

	##number of random samples to use when estimating delta
	delta_rand_samples = 10000;
	##significance threshold for delta
	alpha = 0.01;
	##number of neighbors in the local homogeneity field
	k = 8;

	##load data for domain identification
	data,latitudes,longitudes = load_data(path_to_data, 'sst','lat','lon');
	##normalize time series to zero mean and unit variance
	data = remove_mean(data);
	data = remove_variance(data);

	##convert data to a numpy array where masked values are set to NaN
	data = utils.masked_array_to_numpy(data);
	np.save("test_data",data);
	np.save("lats",latitudes.data);
	np.save("lons",longitudes.data);
	##estimate delta
	print('Estimating delta');
	delta = utils.estimate_delta(data, delta_rand_samples, alpha);
	print('Delta estimate: '+str(delta));
	##step 1. seed identification
	local_homogeity_field, seed_positions = seed_identification(data,latitudes.data,
																longitudes.data,
																delta, k);
	plt.figure();
	plt.imshow(local_homogeity_field-seed_positions);
	plt.savefig('seed_map');
