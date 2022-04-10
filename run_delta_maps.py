import sys
import os
import time
import numpy as np
import netCDF4
import json
import shutil
import argparse
from netCDF4 import Dataset
from Preprocessing.preprocessing import remove_mean
from Preprocessing.preprocessing import remove_variance
from Utils import utils
from SeedIdentification.seed_identification import seed_identification
from DomainIdentification.DomainIdentification import DomainIdentification
from net_inference.net_inference import net_inference
from matplotlib import pyplot as plt


# Authors
# Ilias Fountalis (Foudalisi@hotmail.com) & Fabrizio Falasca (fabrifalasca@gmail.com)
# License: GNU General Public License v3.0

import time

start = time.time()

# After installing new version of numpy we have warnings.
# If you want to remove them simply uncomment these 2 lines.
#from warnings import simplefilter
#simplefilter(action='ignore', category=DeprecationWarning)


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
    nc_fid = Dataset(path_to_data, 'r')
    climate_field = nc_fid.variables[climate_variable][:]
    latitudes = nc_fid.variables[latitude_string][:]
    longitudes = nc_fid.variables[longitude_string][:]

    return climate_field, latitudes, longitudes


def load_config(path_to_config_file):
    """
    loads the configuration file
    """
    with open(path_to_config_file) as f:
        config = json.load(f)
    return config


'''
def plot_domain_map(domains):
    dmap = np.zeros((domains[0].map.shape[0], domains[0].map.shape[1]));
    for d in domains:
        dmap[d.map == 1] = np.random.randint(1, len(domains)+1);
    dmap[dmap == 0] = None
    plt.figure(figsize=[20,20]);
    plt.imshow(np.flipud(dmap));
    plt.savefig('domain_map');
'''


def init_directories(output_dir):

    if(output_dir[-1] != '/'):
        output_dir += '/'

    if(os.path.exists(output_dir)):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    seed_results_dir = output_dir + "seed_identification/"
    domain_results_dir = output_dir + "domain_identification/"
    network_results_dir = output_dir + "network_inference/"

    os.mkdir(seed_results_dir)
    os.mkdir(domain_results_dir)
    os.mkdir(network_results_dir)

    return seed_results_dir, domain_results_dir, network_results_dir


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", "--input", required=True,
                           help="Path to .json config file")

    args = vars(argparser.parse_args())
    # load config file
    config = load_config(args["input"])

    # init output directory and directories to save results
    seed_results_dir, domain_results_dir, network_results_dir = init_directories(
        config["output_directory"])

    # convert this to user input when done
    path_to_data = config["path_to_data"]

    # number of random samples to use when estimating delta
    delta_rand_samples = config["delta_rand_samples"]
    # significance threshold for delta
    alpha = config["alpha"]
    # number of neighbors in the local homogeneity field
    k = config["k"]

    # load data for domain identification
    data, latitudes, longitudes = load_data(path_to_data, config["variable_name"],
                                            config["latitude_name"], config["longitude_name"])

    # normalize time series to zero mean and unit variance
    data = remove_mean(data)
    data = remove_variance(data)

    # convert data to a numpy array where masked values are set to NaN
    data = utils.masked_array_to_numpy(data)

    # estimate delta
    print('Estimating delta')
    delta = utils.estimate_delta(data, delta_rand_samples, alpha);
    print('Delta estimate: '+str(delta))

    # step 1. seed identification
    local_homogeneity_field, seed_positions = seed_identification(data, latitudes.data,
                                                                  longitudes.data,
                                                                  delta, k)
    np.save(seed_results_dir+"local_homogeneity_field", local_homogeneity_field)
    np.save(seed_results_dir+"seed_positions", seed_positions)

    # step 2. domain identification
    print('Domain identification started')
    domain_identifier = DomainIdentification(data, latitudes, longitudes, seed_positions,
                                             k, delta, domain_results_dir)
    domain_identifier.domain_identification()
    domains = domain_identifier.domains
    # domain_identifier.dump_output();

    # From the domains object let's retrieve the domains maps and ids
    # How many domains
    length_dom = len(domains)

    # Let's retrieve the domains' maps and ids
    dom_maps = []
    dom_ids = []
    for i in range(length_dom):
        dom_maps.append(domains[i].map)
        dom_ids.append(domains[i].domain_id)
    dom_maps = np.array(dom_maps)
    dom_ids = np.array(dom_ids)

    print('Done with domain identification')

    # Save them in the output folder
    np.save(domain_results_dir+"domain_maps", dom_maps)
    np.save(domain_results_dir+"domain_ids", dom_ids)

    print('Done')

    # step 3. network inference
    # load data again for network inference
    data, latitudes, longitudes = load_data(path_to_data, config["variable_name"],
                                            config["latitude_name"], config["longitude_name"])

    print('')
    print('Network inference started')
    # Infer the netork
    network_list, strength_list, strength_map = net_inference(
        dom_maps, dom_ids, data, latitudes, longitudes, config["tau_max"], config["q"])
    print('Done')
    print('')
    print('Saving the outputs')
    np.save(network_results_dir+"network_list", network_list)
    np.save(network_results_dir+"strength_list", strength_list)
    np.save(network_results_dir+"strength_map", strength_map)

    end = time.time()
    print('Finished in '+str(round(end - start, 2))+' seconds')
