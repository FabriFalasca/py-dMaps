import sys
import os
import time
import numpy as np
from sklearn.neighbors import NearestNeighbors
from Utils.utils import lagged_correlation
#from matplotlib import pyplot as plt


def seed_identification(data, latitudes, longitudes, delta, k):

	##convert to radians
	latitudes_rad = np.radians(latitudes.copy());
	longitudes_rad = np.radians(longitudes.copy());
	##increase k by 1, since we are using nearest neighbors we always include our self
	k+=1; 

	coords = [];
	indices = [];
	for i,lat in enumerate(latitudes_rad):
		for j,lon in enumerate(longitudes_rad):
			coords.append([lat,lon]);
			indices.append([i,j]);
	coords = np.asarray(coords,dtype=np.float32);
	indices = np.asanyarray(indices,dtype = np.int32);


	nearest_neighbors = NearestNeighbors(n_neighbors=k, algorithm = 'ball_tree',metric='haversine');
	nearest_neighbors.fit(coords);

	##construct local homogeneity field
	print('Calculating local homogeneity field')
	counter = 0;
	local_homogeneity_field = np.zeros((data.shape[1], data.shape[2]));

	for index_i,index_j in indices:
		##get coordinates of grid point
		point_coords = np.expand_dims(coords[counter],0);
		##if the grid cell is masked do not calculate homogeneity field
		if(np.isnan(data[0,index_i,index_j])):
			local_homogeneity_field[index_i, index_j] = np.nan;
		else:
			##coordinates of the k closest points
			_,k_neighborhood_indices = nearest_neighbors.kneighbors(point_coords);
			k_neighborhood_indices = k_neighborhood_indices[0];

			##calculate average correlation at 0 lag between all pairs of indices
			local_homogeneity = 0;
			valid_pairs = 0;
			for i in range(len(k_neighborhood_indices)):
				indices_i = indices[k_neighborhood_indices[i]];
				timeseries_i = data[:, indices_i[0],indices_i[1]];
				if(np.isnan(timeseries_i[0])):
					continue;
				for j in range(i+1,len(k_neighborhood_indices)):
					indices_j = indices[k_neighborhood_indices[j]];
					timeseries_j = data[:, indices_j[0],indices_j[1]];
					if(np.isnan(timeseries_j[0])):
						continue;
					corr = lagged_correlation(timeseries_i, timeseries_j,tau=0,normed=True);
					local_homogeneity += corr;
					valid_pairs += 1;
			##get the average correlations (local homogeneity)
			if(valid_pairs == 0): ##grid cell inland
				local_homogeneity = np.nan;
			else:
				local_homogeneity = local_homogeneity/valid_pairs;

			local_homogeneity_field[index_i, index_j] = local_homogeneity;

		counter += 1;
    
		if(counter%1000 == 0):
			print('Progress: '+str(np.round(counter/len(indices),2)));




	##identify seeds
	print('Seed identification.')
	seed_positions = np.zeros((data.shape[1], data.shape[2]));
	counter = 0;

	for index_i, index_j in indices:
		if(np.isnan(data[0,index_i,index_j])):
			pass;
		else:
			point_coords = np.expand_dims(coords[counter],0);
			_,k_neighborhood_indices = nearest_neighbors.kneighbors(point_coords);
			k_neighborhood_indices = k_neighborhood_indices[0];
			k_neighborhood_indices = k_neighborhood_indices[1:]; ##remove self
        
			##list that holds the homogeneities of the grid cells neighbors
			cell_homogeneities = [];
			for idx in k_neighborhood_indices:
				cell_idx = indices[idx];
				cell_homogeneity = local_homogeneity_field[cell_idx[0],cell_idx[1]];
				cell_homogeneities.append(cell_homogeneity);
			##first entry in the cell_homogeneities is the grid self itself, so easy way to 
			##see if it is a local maximum AND larger than delta
			seed_homogeneity = local_homogeneity_field[index_i, index_j];
			if(seed_homogeneity > max(cell_homogeneities) and seed_homogeneity > delta):
				seed_positions[index_i,index_j] = 1;
		counter += 1;
		if(counter%1000 == 0):
			print('Progress: '+str(np.round(counter/len(indices),2)));


	return local_homogeneity_field, seed_positions;
