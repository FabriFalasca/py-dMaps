import sys
import os
import time
import numpy as np
from sklearn.neighbors import NearestNeighbors
from Utils.utils import lagged_correlation
from DataStructures.GridCell import GridCell
from DataStructures.Domain import Domain
from DataStructures.OverlappingDomain import OverlappingDomain


class DomainIdentification(object):

	def __init__(self, data, latitudes, longitudes, seed_positions, k, delta, output_dir):
		self.data = data;
		self.latitudes = np.radians(latitudes);
		self.longitudes = np.radians(longitudes);
		self.seed_positions = seed_positions;
		self.k = k;
		self.delta = delta;
		self.output_dir = output_dir;
		
		##init. coordinates (latitudes and longitudes) and indices (x,y positions)
		self.coords, self.indices = self._init_coords_and_indices();
		##init. nearest neighbors search
		self.nearest_neighbors = self._init_nearest_neighbors(self.k);
		##init. grid cell dict and grid cell objects and mapping from x,y, coordinate to grid cell id
		self.grid_cells, self.pos_to_grid_cell_id = self._init_grid_cells();
		##init. domains
		self.domains = self.domain_initialization();


	def domain_identification(self):
		"""
		Main function that starts domain identification
		"""
		iteration = 0;
		while(True): ##while we can either expand or merge domains
			##step 1. merge domains
			print('Merging domains');
			merged = self.merge_domains();
			##step 2. expand domains
			expanded, merged = self.expand_domains();

			print('Iteration: '+str(iteration)+", expanded = "+str(expanded)+", merged: "+str(merged));
			iteration += 1;
			if(not merged and not expanded):
				break;



	def merge_domains(self):
		"""
		Somain merging op.
		"""
		merged = False;
		while(True): ##while we can merge domains...
			##get all pairs of overlapping domains
			print('Searching for overlaping domains')
			print('Current domains: '+str(len(self.domains)));

			overlapping_domains = [];
			for i in range(len(self.domains)):
				for j in range(i+1,len(self.domains)):
					if(self.does_overlap(self.domains[i], self.domains[j])):
						overlapping_domain = OverlappingDomain(self.domains[i].domain_id,
																self.domains[j].domain_id);
						overlapping_domain.map = self.domains[i].map + self.domains[j].map;
						overlapping_domains.append(overlapping_domain);

			if(len(overlapping_domains) == 0):
				break; ## no overlapping domains found, stop. 

			##calculate the homogeneity of the overlapping domains.
			for ov_domain in overlapping_domains:
				ov_domain.homogeneity = self.get_homogeneity_overlapping_domain(ov_domain);

			##sort pverlapping domains by homogeneity in deescending order
			overlapping_domains.sort(key=lambda x: x.homogeneity, reverse=True);

			##get the best overlapping domain (max homogeneity)
			best_ov_domain = overlapping_domains[0];

			if(best_ov_domain.homogeneity > self.delta): ##merge
				print('Overlapping domain found, homogeneity: '+str(best_ov_domain.homogeneity));
				##remove the two domains that need to be merged
				idx_a = self.domains.index(Domain(best_ov_domain.id_a));
				domain_a = self.domains.pop(idx_a);
				idx_b = self.domains.index(Domain(best_ov_domain.id_b));
				domain_b = self.domains.pop(idx_b);

				##construct the new domain...
				new_domain = Domain(domain_a.domain_id); ##randomly assign one id
				new_domain.homogeneity = best_ov_domain.homogeneity;

				##add grid cells to the new domain
				for grid_cell_id, grid_cell in domain_a.grid_cells.items():
					new_domain.grid_cells[grid_cell_id] = grid_cell;
				for grid_cell_id, grid_cell in domain_b.grid_cells.items():
					new_domain.grid_cells[grid_cell_id] = grid_cell;
				##update the map
				new_domain.map = np.zeros((domain_a.map.shape[0],domain_a.map.shape[1]));
				new_domain.map[domain_a.map == 1] = 1;
				new_domain.map[domain_b.map == 1] = 1;
				self.domains.append(new_domain);
				merged = True;
			else:
				break;
		return merged;




	def get_homogeneity_overlapping_domain(self, ov_domain):
		"""
		homogeneity for an overlapping domain object
		"""
		indices = np.argwhere(ov_domain.map > 0);
		homogeneity = 0;
		N = len(indices);
		##compute correlation between all pairs of grid cells
		for i in range(N):
			ts_a = self.data[:,indices[i][0],indices[i][1]];
			for j in range(i+1, N):
				ts_b = self.data[:,indices[j][0],indices[j][1]];
				corr = lagged_correlation(ts_a, ts_b,tau=0,normed=True);
				homogeneity += corr;

		##normalize
		n_pairs = N*(N-1)/2.;
		homogeneity = homogeneity/n_pairs;
		return homogeneity;



	def expand_domains(self):
		"""
		Main function to expand domains, it tries to expand all domains
		starting from the highest homogeneity ones to the lowest one.
		At each step it checks if after expansion merging is possible
		"""

		expanded = False;
		start_merging = False;

		while(not start_merging):
			self.domains.sort(key=lambda x: x.homogeneity, reverse=True);

			##for each domain
			for i, domain in enumerate(self.domains):
				domain, expanded = self.expand_domain(domain);
				if(expanded):
					"""
                	Main logic.
                	if the domain is expanded try to see if we have any overlapping domains
                	if we have overlapping domains check to see if their union is larger
                	than delta.
                	if union is larger than delta we need to stop and proceed to the
                	expansion step.
                	else keep expanding other domains.
                	"""
					for j, a_domain in enumerate(self.domains):
						if(i!= j): ##do not check with overlap with yourself
							if(self.does_overlap(domain, a_domain)):
								##these domains overlap, check homogeneity of union
								if(self.get_homogeneity_union(domain,a_domain) > self.delta):
									start_merging = True;
									break;
				if(start_merging):
					break;
			if(not expanded):
				break;

		return expanded, start_merging;




	def get_homogeneity_union(self, domain_a, domain_b):
		"""
		Returns the homogeneity of the union (merging) of two domains
		"""
		##get the overlap map
		overlap_map = domain_a.map+domain_b.map;
		##get the indices of the grid cells of both domains
		indices = np.argwhere(overlap_map > 0);
		homogeneity = 0;

		N = len(indices);

		##compute correlation between all pairs of grid cells
		for i in range(N):
			ts_a = self.data[:,indices[i][0],indices[i][1]];
			for j in range(i+1, N):
				ts_b = self.data[:,indices[j][0],indices[j][1]];
				corr = lagged_correlation(ts_a, ts_b,tau=0,normed=True);
				homogeneity += corr;
		##normalize
		n_pairs = N*(N-1)/2.
		homogeneity = homogeneity/n_pairs;
		return homogeneity;


	def does_overlap(self,domain_a, domain_b):
		##check if two domains overlap
		overlap_map = domain_a.map+domain_b.map;
		overlaps = overlap_map[overlap_map>1];
		if(len(overlaps) > 0):
			return True;
		else:
			return False;

	def expand_domain(self, domain):
		"""
		Performs one step of domain expansion.
		Returns:
			domain: The expanded domain (if domain could be expanded)
			is_expanded: (boolean) True if domain is expanded False if not
		"""
		is_expanded = False;

		##generate a set that will have all grid cell ids that are
		##already inside the domain
		domain_grid_cells = set();
		for grid_cell_id, grid_cell in domain.grid_cells.items():
			domain_grid_cells.add(grid_cell);

		##set that will hold all candidate grid cells.
		candidate_grid_cells = set(); ##list with GridCell objects

		##get candidate grid cells for expansion
		for grid_cell_id, grid_cell in domain.grid_cells.items():
			##for each grid cell get its k=8 nearest neighbors
			coords = np.expand_dims(np.asarray([grid_cell.lat,grid_cell.lon]),0);
			_,k_neighborhood_indices = self.nearest_neighbors.kneighbors(coords,n_neighbors=8+1);
			k_neighborhood_indices = k_neighborhood_indices[0];
			##remoove self
			k_neighborhood_indices = k_neighborhood_indices[1:];

			for k_index in k_neighborhood_indices:
				neigh_indices = self.indices[k_index]; ##x,y index of neighboring grid cell 
				##get the corresponding grid cell
				try:
					##neighboring grid cell id
					neigh_grid_cell_id = self.pos_to_grid_cell_id[str(neigh_indices[0])+","+str(neigh_indices[1])];
					candidate_grid_cell = self.grid_cells[neigh_grid_cell_id];
					if(candidate_grid_cell not in domain_grid_cells): ##grid cell not in domain = grid cell in border
						candidate_grid_cells.add(candidate_grid_cell);
				except KeyError:
					pass; ##grid cell in land

		##compute average  correlation between each candidate
		##grid cell and grid cells inside the domain

		##check number 1. do we have any candidate grid cells
		if(len(candidate_grid_cells) == 0):
			return domain, is_expanded;

		for candidate_grid_cell in candidate_grid_cells:
			candidate_grid_cell_ts = self.data[:,candidate_grid_cell.index[0],candidate_grid_cell.index[1]];
			score = 0.;

			for domain_grid_cell in domain_grid_cells:
				domain_grid_cell_ts = self.data[:,domain_grid_cell.index[0],domain_grid_cell.index[1]];
				correlation = lagged_correlation(candidate_grid_cell_ts, domain_grid_cell_ts,
													tau=0,normed=True);
				score += correlation;
			score = score/len(domain_grid_cells);
			candidate_grid_cell.score = score;


		##sort grid cells by score in descending order.
		candidate_grid_cells = list(candidate_grid_cells);
		candidate_grid_cells.sort(key=lambda x: x.score, reverse=True);

		##get the best candidate
		best_candidate = candidate_grid_cells[0];
		if(best_candidate.score > self.delta): ##then add the grid cell
			##update domain homogeneity
			domain_pairs = len(domain.grid_cells)*(len(domain.grid_cells)-1)/2.
			current_homogeneity = domain.homogeneity*domain_pairs;
			bc_score = best_candidate.score*len(domain.grid_cells);
			domain.homogeneity = (bc_score+current_homogeneity)/(len(domain.grid_cells)*(len(domain.grid_cells)+1)/2.);
			domain.grid_cells[best_candidate.id] = best_candidate;
			##update domain map
			domain.map[best_candidate.index[0], best_candidate.index[1]] = 1;
			is_expanded = True;
			return domain, is_expanded;
		else: ##this domain can not be expanded
			return domain, is_expanded;


	def domain_initialization(self):
		##get the locations of the seeds
		seed_locations = np.where(self.seed_positions == 1);
		domains = [];

		for i in range(len(seed_locations[0])): ## 1 domain per seed location
			domain = Domain(i); ##init empty domain
			domain.map = np.zeros((self.data.shape[1], self.data.shape[2]),
									dtype = np.int32); ##init domain map;

			##add to the domain map the seed location
			init_index = np.asarray([seed_locations[0][i],seed_locations[1][i]])
			domain.map[init_index[0], init_index[1]] = 1;

			##get the corresponding lat/lon coordinates
			coord_index = np.where( (self.indices[:,0] == init_index[0]) &
									(self.indices[:,1] == init_index[1])
								  )[0];
			init_coords = self.coords[coord_index]; ##corresponding lat,lon of init_index

			##and no go and get the k nearest neighbors and add them as well
			_,k_neighborhood_indices = self.nearest_neighbors.kneighbors(init_coords);
			k_neighborhood_indices = k_neighborhood_indices[0];

			for k_index in k_neighborhood_indices:
				neigh_indices = self.indices[k_index];
				##get the corresponding grid cell
				try:
					grid_cell_id = self.pos_to_grid_cell_id[str(neigh_indices[0])+","+str(neigh_indices[1])];
					grid_cell = self.grid_cells[grid_cell_id];
					##add it to the domain
					domain.grid_cells[grid_cell_id] = grid_cell;
					##update domain map
					domain.map[neigh_indices[0], neigh_indices[1]] = 1;
				except KeyError:
					pass; ##grid cell in land

			##final step, initialize the local homogeneity of the domain
			local_homogeneity = 0;
			valid_pairs = 0;
			for i in range(len(k_neighborhood_indices)):
				indices_i = self.indices[k_neighborhood_indices[i]];
				timeseries_i = self.data[:, indices_i[0],indices_i[1]];
				if(np.isnan(timeseries_i[0])): ##in land...
					continue;
				for j in range(i+1,len(k_neighborhood_indices)):
					indices_j = self.indices[k_neighborhood_indices[j]];
					timeseries_j = self.data[:, indices_j[0],indices_j[1]];
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
			domain.homogeneity = local_homogeneity;
			domains.append(domain);

		return domains;

	def _init_grid_cells(self):
		grid_cells = {}; ##key grid cell id, value: GridCell
		##maps a grid cell coordinate (x,y location) to a grid cell id
		pos_to_grid_cell_id = {}; ##key str ("x,y") position, value: grid cell id.
	
		grid_cell_id = 0; ##initial grid cell id
		for i in range(self.latitudes.shape[0]):
			for j in range(self.longitudes.shape[0]):
				if(not np.isnan(self.data[0,i,j])): ##grid cell not in land
					grid_cell = GridCell(grid_cell_id);
					grid_cell.index[0] = i;
					grid_cell.index[1] = j;
					grid_cell.lat = self.latitudes[i];
					grid_cell.lon = self.longitudes[j];
					grid_cells[grid_cell_id] = grid_cell;

					pos_to_grid_cell_id[str(i)+","+str(j)] = grid_cell_id;
					grid_cell_id += 1;
		return grid_cells,pos_to_grid_cell_id;

	
	def _init_nearest_neighbors(self, k):
		nearest_neighbors = NearestNeighbors(n_neighbors=k+1, algorithm = 'ball_tree',metric='haversine');
		nearest_neighbors.fit(self.coords);
		return nearest_neighbors;
	
	def _init_coords_and_indices(self):

		coords = []; ##array with tuple latitude, longitude
		indices = []; ##array with index (i,j) of grid cell

		for i, lat in enumerate(self.latitudes):
			for j, lon in enumerate(self.longitudes):
				coords.append([lat,lon]);
				indices.append([i,j]);
		coords = np.asarray(coords,dtype=np.float32);
		indices = np.asanyarray(indices,dtype = np.int32);
		return coords, indices
	

	def dump_output(self):
		"""
		For each domain it creates:
			(a) a numpy array with a mask set to 1 if grid cells exist in the x,y location of the map
			(b) a text file with raw time series in the form x,y,ts[0],...,[ts[T-1]
			(c) a file with domain information such as homogeneity, number of grid cells etc
		"""
		for domain in self.domains:
			np.save(self.output_dir+"map_"+str(domain.domain_id),domain.map);
			
			with open(self.output_dir+"info_"+str(domain.domain_id)+".txt",'w') as f:
				f.write("Number_of_grid_cells:"+str(len(domain.grid_cells))+"\r\n");
				f.write("Homogeneity:"+str(domain.homogeneity)+"\r\n");
				f.write("Delta_est:"+str(self.delta)+"\r\n");
				f.write("K_was_set_to:"+str(self.k)+"\r\n");



