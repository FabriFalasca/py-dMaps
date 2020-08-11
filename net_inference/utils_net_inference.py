import sys
import os
import time
import numpy as np
import scipy.stats


def masked_array_to_numpy(data):
    return np.ma.filled(data.astype(np.float32), np.nan);

def get_nonmask_indices(data):
    return np.argwhere(~np.isnan(np.sum(data,axis=0)));

# Covariance between two time series x(t) and y(t) at lag tau
def lagged_covariance(x,y,tau,normed=False):

    assert len(x) == len(y);
    #length of time series
    T = len(x);

    ##assert that lag can nit be greater than T
    assert tau < T;

    ##if time series are not normalized, reduce them to zero mean and unit variance
    if(not normed):
        x = x-np.mean(x);
        y = y-np.mean(y);


    if(tau == 0):
        return np.dot(x,y)/T;
    if(tau > 0):
        x = x[0:T-tau];
        y = y[tau:T]
        return np.dot(x,y)/T;

    if(tau < 0):
        tau = np.abs(tau)
        y = y[0:T-tau];
        x = x[tau:T]
        return np.dot(x,y)/T;

# Get the covariances between timeseries ts1 and ts2 for a range of lags
def get_covariances(ts1,ts2,tau_range,normed=False):
    assert len(ts1) == len(ts2);
    covariances = [];
    for tau in range(-tau_range,tau_range+1):
        covariances.append(lagged_covariance(ts1,ts2,tau,normed));
    return covariances;

# Correlation between two time series x(t) and y(t) at lag tau
def lagged_correlation(x,y,tau,normed=False):

    assert len(x) == len(y);
    #length of time series
    T = len(x);

    ##assert that lag can nit be greater than T
    assert tau < T;


    ##if time series are not normalized, reduce them to zero mean and unit variance
    if(not normed):
        x = (x-np.mean(x))/np.std(x);
        y = (y-np.mean(y))/np.std(y);


    if(tau == 0):
        return np.dot(x,y)/T;
    if(tau > 0):
        #corr = 0; Changed by Fab
        x = x[0:T-tau];
        y = y[tau:T]
        return np.dot(x,y)/T;

    if(tau < 0):
        tau = np.abs(tau)
        y = y[0:T-tau];
        x = x[tau:T]
        return np.dot(x,y)/T;

# Autocorrelation of timeseries ts at lag tau
def lagged_autocorrelation(ts,tau,normed=False):
    return lagged_correlation(ts,ts,tau,normed);

# Get the correlations between timeseries ts1 and ts2 for a range of lags
def get_correlogram(ts1,ts2,tau_range,normed=False):
    assert len(ts1) == len(ts2);
    correlogram = [];
    for tau in range(-tau_range,tau_range+1):
        correlogram.append(lagged_correlation(ts1,ts2,tau,normed));
    return correlogram;

# Bartlett variance at lat tau
def bartlett_variance_tau(ts1,ts2,tau,normed=False):

    assert len(ts1) == len(ts2);
    T = len(ts1);
    ##get the correlogram of the two time series
    correlogram_ts1 = get_correlogram(ts1,ts1,T-1,normed);
    correlogram_ts2 = get_correlogram(ts2,ts2,T-1,normed);

    var = np.sum(np.multiply(correlogram_ts1,correlogram_ts2))/(T-tau);
    # set to zero (small) negative numbers
    if(var <= 0):
        var = np.random.uniform(0, 0.000001)
    return var;

# Function for the edge inference.
def edge_inference(tau_max,corr_matrix,cov_matrix,bartlett_std_matrix,domains_pair):

    # input:
    #        (a) tau_max: it defines the lag range R \in [-tau_max,tau_max]
    #        (b) corr_matrix: correlograms for each pair of timeseries
    #        (c) cov_matrix: covariances for range of lags for each pair of timeseries
    #        (d) bartlett_std_matrix: bartlett standard deviations for each pair of timeseries
    #        (e) all pairs of domains (same order of correlograms)

    # output list of edges.
    #        Each edge have the following format
    #        [domain a, domain b, tau_min, tau_max, tau*, r*, weight]
    #        where

    #       r* is the max significant correlation
    #       tau* is the lag correspondent to r*

    #       [tau_min,tau_max] defines the range of lags associated with significant correlations
    #       in the interval [r*-bartlett std(r*),r*+bartlett std(r*)]
    #
    #       if the range [tau_min, tau_max] include zero: domain a <---> domain b
    #       if tau_min > 0:                               domain a  ---> domain b
    #       if tau_max < 0:                               domain a <---  domain b

    #       weight is the link weight: covariance at tau*


    # Define an array with the lags
    lags = np.arange(-tau_max,tau_max+1,1)

    # initialize the array that will hold the edges
    edges = []

    for i in range(len(corr_matrix)):

        # If the correlogram is all zeros -> no signifinant correlations -> skip it
        if np.sum(np.abs(corr_matrix[i])) == 0:

            continue

        else:
            # Correlograms of link domains_pair[i]
            correlogram = corr_matrix[i]
            # Covariances of link domains_pair[i]
            covariance = cov_matrix[i]
            # Bartletts standard deviations of the correlogram relative to domains_pair[i]
            bartlett_sigmas = bartlett_std_matrix[i]

            # Compute the max in absolute value
            max_corr = max(correlogram,key=abs)
            # Find the correspondent bartlett std
            std_max = bartlett_sigmas[correlogram==max_corr][0]
            # Keep only the values that are inside 1 sigma of the max
            constraints = np.logical_and(correlogram<=max_corr+std_max, correlogram>=max_corr-std_max)
            # Find the lags associated with the correlations respecting correlogram
            lags_connection = lags[constraints]
            # Min and max lag
            lag_min = lags_connection[0]
            lag_max = lags_connection[-1]
            # Lag tau* correspondent to r*
            lag_star = lags[correlogram==max_corr][0]
            # Weight: covariance at Lag tau*
            weight = covariance[correlogram==max_corr][0]

            edge = [domains_pair[i][0],domains_pair[i][1],lag_min,lag_max,lag_star,max_corr,weight]

        edges.append(edge)

    return edges

# Network inference using FDR
def net_inference_FDR(path,ids,tau_max,q):
    # Inputs:
    # path to time series
    # tau_max: it defines the lag range R \in [-tau_max,tau_max]
    # q: false discovery rate
    # list of ids of domains (i.e., ids = ['0', '20', '22', '35', ...])

    # Outputs:
    # list with all connections
    # list with all strengths

    # import the time series
    ts = np.loadtxt(path)

    # Number of signals
    N = len(ts)

    # the time series imported have zero mean but NOT 1 sigma
    # We standardize the ts for the correlation analysis:

    # mean(normed_ts[i]) = 0 && std(normed_ts[i]) = 1 for all i
    normed_ts = (ts.T/np.std(ts,axis=1)).T

    ## Compute the correlograms for each unique pair
    correlograms = []
    for i in range(N):
        ts1 = normed_ts[i]
        for j in range(i+1,N):
            ts2 = normed_ts[j]
            correlograms.append(get_correlogram(ts1,ts2,tau_max,normed=True))

    correlograms = np.asarray(correlograms)

    ## Compute the lag covariances for each unique pair
    covariances = []
    for i in range(N):
        ts1 = ts[i]
        for j in range(i+1,N):
            ts2 = ts[j]
            covariances.append(get_covariances(ts1,ts2,tau_max,normed=True))

    covariances = np.asarray(covariances)

    ## We want the Bartlett's variance for each unique pair of time series
    bartlett = np.zeros([int(N*(N-1)/2),int(2*tau_max+1)])
    k1 = 0
    for i in range(N):
        ts1 = normed_ts[i]
        for j in range(i+1,N):
            ts2 = normed_ts[j]
            k2 = 0
            for tau in np.arange(-tau_max,tau_max+1,1):
                bartlett[k1,k2] = bartlett_variance_tau(ts1,ts2,tau,normed=True)
                k2 += 1
            k1 += 1

    ## Compute the z score for every correlation

    # By taking the absolute value of the correlations is like doing a 2 - tailed t-test
    abs_correlograms = np.abs(correlograms)
    bartlett_std = np.sqrt(bartlett)
    z_scores = abs_correlograms/bartlett_std

    # Compute all p-values for each connection
    p_vals = np.zeros([int(N*(N-1)/2),int(2*tau_max+1)])

    for i in range(np.shape(z_scores)[0]):
        for j in range(np.shape(z_scores)[1]):
            ##one sided t-test
            p_vals[i,j] = 1 - scipy.stats.norm(0,1).cdf(z_scores[i,j]);

    # Compute p_min

    # Step (a)
    # flatten all p_values and sort them in descending order
    sorted_p_vals = sorted(p_vals.flatten())
    sorted_p_vals = np.asarray(sorted_p_vals)
    # How many p-values
    n_p_vals = len(sorted_p_vals)
    # Step (b)
    # Define the line for computing the p_min
    fdr_ratio = (q/n_p_vals)*np.arange(1,n_p_vals+1)

    ## Get p_min
    # (a) take the difference between fdr_ratio[i] and sorted_p_vals[i]
    difference = fdr_ratio-sorted_p_vals
    # (b) take only the positive entries
    positive_vals = sorted_p_vals[difference>0]
    # (c) p_min is the last entry
    p_min = positive_vals[-1]

    ## pairs is an array with all possible pairs
    # It has the same order of the arrays correlations or covariances

    # E.g., correlations[0] is the correlogram between the pair of domains in pairs[0]

    pairs = []
    for i in range(N):
        id1 = ids[i]
        for j in range(i+1,N):
            id2 = ids[j]
            pairs.append([id1,id2])

    pairs = np.array(pairs)

    ## We set to zero all the correlations with p value > p_min
    correlograms[p_vals>p_min] = 0

    ## Compute the network
    network = edge_inference(tau_max,correlograms,covariances,bartlett_std,pairs)

    # Compute the strengths
    indices_net = np.array(network)[:,0:2]

    strength_list = []

    for i in range(len(ids)):

        # id considered
        id_domain = ids[i]

        if len(np.where(indices_net==ids[i])[0]) == 0:

            strength_list.append([id_domain,0])

        else:

            # Positions of the index in the net
            pos = np.where(indices_net==ids[i])[0]

            # Initialize a sublist where we hold the weights associated to each connections
            sublists = []

            for j in range(len(pos)):

                sublists.append(network[pos[j]][6]) # the weight is entry 6 in the network array

            strength_list.append([id_domain,np.sum(np.abs(sublists))])

    return network, strength_list
