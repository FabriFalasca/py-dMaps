# -*- coding: utf-8 -*-
# importing libraries

import numpy as np
import scipy.io
import numpy.ma as ma
from scipy import signal, spatial
from scipy.spatial.distance import pdist,cdist, squareform
from collections import Counter

# import utils
import utils
# Network libraries
import networkx as nx
import infomap
import preprocessing

# Fabri Falasca
# fabrifalasca@gmail.com

'''
Given time series embedded in a spatiotemporal grid,
this code find clusters/modes based on their correlation
networks. The clustering problem is cast as a community detection
problem.
The community detection algorithm used is infomap
'''

# Function for Pearson correlation
def pearson_corr(x, y):
    # Time series are supposed to be normalized to
    # (a) zero mean and (b) unit variance
    # Length of time series
    T = len(x)
    return np.dot(x,y)/T

# Function to compute the correlation matrix using the Pearson correlation
# defined above
def compute_corr_m(xs, ys):
    return cdist(xs, ys, metric=pearson_corr)

# Function to compute the Adjacency matrix from the correlation matrix
def compute_adj_m(x,k):
    # inputs: - x is a correlation (distance) matrix
    #         - tau is a threshold
    adj_m = x.copy()
    # Step (1) set the diagonal of the correlation matrix to zero
    np.fill_diagonal(adj_m,0)
    # Step (2) for a given number i, if i >= tau ---> 1, if i < tau ---> 0
    adj_m[adj_m>=k]=1
    adj_m[adj_m<k]=0

    return adj_m

# Function to compute the adjacency matrix weighted
def compute_adj_m_weighted(x):
    # inputs: - x is a correlation (distance) matrix
    #         - tau is a threshold
    adj_m = x.copy()
    # Step (1) set the diagonal of the correlation matrix to zero
    np.fill_diagonal(adj_m,0)
    # Step (2)
    # if links are negative assign zero. If not leave the weight
    #adj_m[adj_m>=k]=1
    adj_m[adj_m<0]=0

    return adj_m

'''
Main function

Inputs
(a) path to netcdf file
(b) name of the climate variable (e.g., 'sst')
(c) name for the "longitude" variable (e.g., 'lon')
(d) name for the "latitude" variable (e.g., 'lat')
(e) rand_samples: # of pair of time series for the tau inference
(f) alpha significance level for the tau inference

'''

def infomap_communities(path,climate_variable,lon_variable,lat_variable,
                            rand_samples, alpha):

    # Import data
    data = utils.load_data(path, climate_variable)
    lon = utils.load_data(path, lon_variable)
    lat = utils.load_data(path, lat_variable)

    ##normalize time series to zero mean and unit variance
    data = preprocessing.remove_mean(data)
    data = preprocessing.remove_variance(data)

    # Dimension of the field
    dimX = np.shape(data)[2] # Points in X axis (long)
    dimY = np.shape(data)[1] # Points in Y axis (lat)
    dimT = np.shape(data)[0] # Points in time

    # Check: if a grid point is masked with nan at time t it should be masked at all t
    for i in range(dimY):
        for j in range(dimX):
            if np.isnan(np.sum(data[:,i,j])):
                data[:,i,j] = np.nan
            else:
                data[:,i,j] = data[:,i,j]

    # From mask array to np.array
    data = utils.masked_array_to_numpy(data)
    lon = utils.masked_array_to_numpy(lon)
    lat = utils.masked_array_to_numpy(lat)

    # To remove a trend we high pass filter the time series
    #data = high_pass_field(data)

    # From a [dimT,dimY,dimX] array
    # to a list of (dimY x dimX) spectra of length dimV
    flat_data_masked = data.reshape(dimT,dimY*dimX).transpose()

    # Consider only the points that are not masked
    flat_data = flat_data_masked[~np.isnan(np.sum(flat_data_masked,axis=1))]

    # Compute tau
    # function "estimate K" is in utils.py
    print('Starting K inference')
    print('Random sample set to '+str(rand_samples))
    print('Significance for tau inference set to alpha = '+str(alpha))
    k, corrs = utils.estimate_k(data, rand_samples, alpha)
    print('K = '+str(k))

    # Now that we have tau let's compute the adjacency matrix
    # Compute the correlation matrix (this is the longest thing to compute and store)
    print('Computing the correlation matrix')
    corr_matrix = compute_corr_m(flat_data, flat_data)

    # Compute the adjacency matrix A
    # The network is encoded in A
    print('Computing the adjacency matrix')
    adj_matrix = compute_adj_m(corr_matrix,k)
    #adj_matrix = compute_adj_m_weighted(corr_matrix)

    # From an Adjacency matrix to a graph object
    graph = nx.from_numpy_matrix(adj_matrix)

    # Infomap
    print('Community detection via Infomap')

    im = infomap.Infomap("--two-level --verbose --silent")

    # Add linnks
    for e in graph.edges():
        im.addLink(*e)

    # run infomap
    im.run()

    print(f"Found {im.num_top_modules} communities")

    # List of nodes ids
    nodes = []
    # List of respective community
    communities = []

    for node, module in im.modules:

        nodes.append(node)
        communities.append(module)

    partition = np.transpose([nodes,communities])

    # How many communities?
    n_com = np.max(partition[:,1])

    # Time series should be weighted by the cosine
    # of their latitude.
    print('Weighting time series')
    # Transform latitudes in radians
    latitudes = np.radians(lat);
    # Assign a weight to each latitude phi
    lat_weights = np.cos(latitudes).reshape(len(latitudes),1)
    # Define the weighted domain
    data = data * lat_weights
    # Let's define again the flat_data array
    flat_data_mask_weighted = data.reshape(dimT,dimY*dimX).transpose()
    flat_data_weighted = flat_data_mask_weighted[~np.isnan(np.sum(flat_data_mask_weighted,axis=1))]
    print('Computing signals')
    # In modes_indices we store the indices of nodes belonging to each community
    #  e.g., modes_indices[0] will contain all nodes in community 0 etc.
    modes_indices = []
    for i in np.arange(1,n_com+1):
        modes_indices.append(partition[partition[:,1]==i][:,0])
    modes_indices = np.array(modes_indices)

    # Now compute the signals of each community as the cumulative anomaly inside
    signals = []
    for i in range(len(modes_indices)):
        # Extract the mode
        extract_mode = np.array(flat_data_weighted)[modes_indices[i]]
        # Compute the signal
        signal = np.sum(extract_mode,axis=0)
        # Store the result
        signals.append(signal)

    signals = np.array(signals)


    print('Embed communities in the map')
    
    #########################################################
    ######### This part as to be re-written in a smarter way. 
    #########################################################

    # Back into the map

    # Find the mapping
    # mapping[i] will tell you the mapping from
    # flat_sst[i] ---> k, where k is such that flat_sst_masked[k] == flat_sst[i]
    mapping = []

    # From an array to a list
    list_flat_data_masked = flat_data_masked.tolist()
    list_flat_data = flat_data.tolist()

    for i in range(len(flat_data)):
        mapping.append(list_flat_data_masked.index(list_flat_data[i]))
    mapping = np.array(mapping)

    # Here we store the domains maps (clusters embedded in the grid)
    domains = []

    for i in range(len(modes_indices)):

        # Let's create a temporary spatial grid
        data_grid  = data[0,:,:].copy()
        # Let's flatten it
        flat_data_grid = data_grid.flatten()
        # If is not a nan is a zero
        flat_data_grid[~np.isnan(flat_data_grid)] = 0

        for j in range(len(modes_indices[i])):

            # set to 1 the pixel belonging to a domain
            flat_data_grid[mapping[modes_indices[i][j]]] = 1

        domains.append(flat_data_grid)

    domains = np.array(domains)

    gridded_domains = domains.reshape(len(domains),dimY,dimX)

    # define a single domain map

    domain_map = np.zeros((gridded_domains.shape[1], gridded_domains.shape[2]))
    i = 1
    for d in range(gridded_domains.shape[0]):
        domain_map[gridded_domains[d,:,:] == 1] = i;
        domain_map[np.isnan(gridded_domains[d,:,:])] = np.nan;
        i += 1


    return gridded_domains, signals, domain_map
