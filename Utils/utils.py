import sys
import os
import time
import numpy as np
import scipy.stats




def masked_array_to_numpy(data):
    return np.ma.filled(data.astype(np.float32), np.nan);

def get_nonmask_indices(data):
    return np.argwhere(~np.isnan(np.sum(data,axis=0)));

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


def lagged_autocorrelation(ts,tau,normed=False):
    return lagged_correlation(ts,ts,tau,normed);

def get_correlogram(ts1,ts2,tau_range,normed=False):
    assert len(ts1) == len(ts2);
    correlogram = [];
    for tau in range(-tau_range,tau_range+1):
        correlogram.append(lagged_correlation(ts1,ts2,tau,normed));
    return correlogram;

def bartlett_variance(ts1,ts2,normed=False):

    assert len(ts1) == len(ts2);
    T = len(ts1);
    ##get the correlogram of the two time series
    correlogram_ts1 = get_correlogram(ts1,ts1,T-1,normed);
    correlogram_ts2 = get_correlogram(ts2,ts2,T-1,normed);

    var = np.sum(np.multiply(correlogram_ts1,correlogram_ts2))/T;
    # set to zero (small) negative numbers
    if(var <= 0):
        var = np.random.uniform(0, 0.000001)
    return var;


def estimate_delta(data, rand_samples, alpha):

    ##get all non masked indices
    indices = get_nonmask_indices(data);
    z_scores = [];
    bartlett_vars = [];
    corrs = [];
    for i in range(rand_samples):
        ##sample two time series at random
        idx1 = np.random.randint(0, indices.shape[0],1);
        idx2 = np.random.randint(0, indices.shape[0],1);
        ts1 = data[:, indices[idx1,0], indices[idx1,1]].squeeze();
        ts2 = data[:, indices[idx2,0], indices[idx2,1]].squeeze();

        ##compute variance of time series using bartlett's formula
        bartlett_var = bartlett_variance(ts1,ts2,normed=True);
        bartlett_vars.append(bartlett_var)
        ##compute correlation at zero lag
        corr = lagged_correlation(ts1,ts2,0,normed=True);
        corrs.append(corr);
        ##normalize dot product
        z_score = corr/np.sqrt(bartlett_var);
        z_scores.append(z_score);

        #if(i%100 == 0):
        #    print('Progress: '+str(np.round(i/rand_samples,2)));

    ##compute p-values of z-scores
    p_vals = [];
    for i in range(len(z_scores)):
        ##one sided t-test
        p_val = 1 - scipy.stats.norm(0,1).cdf(z_scores[i]);
        p_vals.append(p_val);
    p_vals = np.asarray(p_vals)
    corrs = np.asarray(corrs);
    ##compute average of all significant correlations
    significant_correlations = corrs[p_vals<alpha];
    print('Number of significant correlations: '+str(len(significant_correlations)))
    delta = np.mean(significant_correlations);

    return delta;
