import sys
import os
import time
import numpy as np
import scipy.stats
import numpy as np
from Utils.utils_net_inference import *
import scipy.stats
import numpy.ma as ma


# Function to infer a network between domains using FDR
def net_inference(dom_maps, dom_ids, data, latitudes, longitudes, tau_max, q):

    # Transform data from "masked array" to "numpy array"
    data = np.ma.filled(data.astype(float), np.nan)
    latitudes = np.ma.filled(latitudes.astype(float), np.nan)
    longitudes = np.ma.filled(longitudes.astype(float), np.nan)

    # Transform latitudes in radians
    latitudes = np.radians(latitudes);
    # Assign a weight to each latitude phi
    lat_weights = np.cos(latitudes).reshape(len(latitudes),1)
    # Define the weighted domain
    weighted_domains = dom_maps*lat_weights

    # Compute all domains signals as the cumulative anomaly inside

    signals = []
    for i in range(len(weighted_domains)):
        signals.append(cumulative_anomaly(data,weighted_domains[i]))
        # signals.append(average_anomaly(data,weighted_domains[i])) # in case you want the average anomaly
    signals = np.array(signals)

    # Compute the network
    network_list, strength_list = net_inference_FDR(signals,dom_ids,tau_max,q)

    strength_list_to_print = strength_list.copy()
    # Compute the strength map
    strength_list = np.array(strength_list)
    # Consider only the strength and not the domains ids
    strengths = strength_list[:,1]

    # Compute the domain' strength map
    dom_strength = []
    for i in range(len(strengths)):
        dom_strength.append(dom_maps[i]*strengths[i])
    dom_strength = np.array(dom_strength)

    # We define the strength map by computing the max strength at each grid
    # point i (domains can overlap)
    strength_map = np.max(dom_strength,axis = 0)

    return network_list, strength_list_to_print, strength_map
