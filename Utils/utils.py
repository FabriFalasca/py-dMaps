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


# Additional functions for domains merging: sometimes useful after running
# the main domain identification algorithm.
# Written by Lucile Ricard (lucile.ricard@epfl.ch)

def sort_all (dom_ids, dom_maps, dom_strength) :
    sort_inds = np.argsort(dom_strength[:,1])
    #Descending order
    sort_inds = np.flip(sort_inds)
    #Sort
    sort_dom_strength = dom_strength[sort_inds]
    sort_dom_ids = dom_ids[sort_inds]
    sort_dom_maps = dom_maps[sort_inds]
    return sort_inds, sort_dom_ids, sort_dom_maps, sort_dom_strength

def sort_domains (dom_ids, dom_maps) :
    d_sizes = []
    for d in dom_maps :
        size = list(d.flatten()).count(1)
        d_sizes.append(size)
    d_sizes = np.array(d_sizes)
    #Sort the domain maps by domain sizes
    inds = np.argsort(d_sizes)
    #In descending order
    inds = np.flip(inds)
    #Sort
    sort_d_sizes = d_sizes[inds]
    sort_dom_ids = dom_ids[inds]
    sort_dom_maps = dom_maps[inds]
    return inds, sort_d_sizes, sort_dom_ids, sort_dom_maps

def test(d1,d2):
    ind1 = np.argwhere(d1 == 1)
    ind2 = np.argwhere(d2 == 1)
    N1 = ind1.shape[0]
    N2 = ind2.shape[0]
    compteur = 0
    for indices in (ind1) :
        i1 = indices[0]
        i2 = indices[1]
        location = np.argwhere(i1 == ind2[:,0])
        if i2 in ind2[location,1] :
            compteur += 1
    if (compteur >= 0.7*N1) & (compteur >= 0.7*N2):
        print('At least 70% of d1 grid cells are d2 grid cells')
        return True
    else :
        return False

def compute_liste (dom_maps) :
    liste = []
    N, dimx, dimy = dom_maps.shape
    for i in range(N) :
        # print('i=', i)
        dmap1 = dom_maps [i,]
        #Read other domain maps still in descending order
        for j in range (i+1,N) :
            # print('j =', j)
            dmap2 = dom_maps[j]
            if test (dmap1,dmap2) == True :
                print('Indices of domain maps to merge = %s and %s' %(i,j))
                liste.append((i,j))
    liste = np.array(liste)
    print(liste)
    return liste

def to_merge (liste, dom_maps, k):
    new_dom_maps = np.copy(dom_maps)
    try :
        i,j = liste [0][0], liste[0][1]
    except :
        print('Empty list. No domains to merge anymore.')
        k +=1
        return new_dom_maps, k
    dmap1 = new_dom_maps[i,]
    dmap2 = new_dom_maps[j,]
    #Take all the 1 from the tow domain maps to merge
    dmap_merged = np.maximum(dmap1, dmap2)
    #Replace the strongest domain map by the merged domain map
    new_dom_maps[i,] = dmap_merged
    #Suppress the other domain map
    new_dom_maps = np.delete(new_dom_maps,j,axis = 0) #length = length -1
    return new_dom_maps, k
