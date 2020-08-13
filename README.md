py-dMAPS

*** THIS README IS UNDER CONSTRUCTION ***

# Contents

(i) Netcdf files 

(ii) Preprocessing

(iii) how to run py-dMAPS

(iv) Outputs

(v) Publications

# (i) Netcdf files 

The code accept netcdf files https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_introduction.html .
Each time series x_i(t) in the dataset is assumed to be already preprocessed (i.e., each x_i(t) should be stationary). In our case we worked often with monthly anomalies: the preprocessing involve (a) trend removal and (b) removal of the seasonal cycle.

To see the structure of a netcdf file "file.nc", open a terminal and type: 

ncdump -h file.nc

To visualize netcdf files please consider the software NCVIEW by David Pierce (freely available at http://meteora.ucsd.edu/~pierce/ncview_home_page.html). Another useful software is Panoply by NASA (freely available at https://www.giss.nasa.gov/tools/panoply/).

# (ii) PREPROCESSING 

To preprocess a netcdf dataset we suggest to use the CDO package from MPI (https://code.mpimet.mpg.de/projects/cdo/).

Below we show how to preprocess a spatiotemporal dataset saved as monthly averages.
Consider a spatiotemporal climate field such as the COBEv2 reanalysis (https://psl.noaa.gov/data/gridded/data.cobe2.html) with temporal and spatial resolution of 1 month and 1 by 1 degree (i.e., 360x180 time series) respectively. The temporal range goes from January 1850 to December 2018.

Possible preprocessing steps COBEv2.nc using CDO:

* Remap the dataset to the resolution you want (i.e, 2 by 2 degree)

cdo -L -remapbil,r180x90 input.nc output.nc

* Select a period (e.g., from January 1950 to December 2015)

cdo selyear,1950/2010 input.nc output.nc

* Remove the seasonal cycle

cdo -L -ymonsub input.nc -ymonmean input.nc output.nc 

* Remove a linear trend

cdo detrend input.nc output.nc

* Focus only on the latitudinal range from 60 degree South to 60 degree North

cdo sellonlatbox,0,360,-60,60 input.nc output.nc

# (iii) HOW TO RUN py-dMAPS 

* python3 run_delta_maps.py -i configs/sample_config.json

Inputs in the configs/sample_config.json

(a) path_to_data: path to the spatiotemporal dataset. Each time series x_i(t) in the dataset is assumed to be already preprocessed (i.e., each x_i(t) should be stationary)

(b) output_directory: name of the output directory that will be created

(c) variable name: name of the variable of interest in the netcdf file

(d) latitude: name of the variable containing latitudes (i.e., lat in our case, but other dataset could have "latitude")

(e) longitude: name of the variable containing longitude (i.e., lon in our case, but other dataset could have "longitude")

(f) delta_rand_samples: random sample of pairs of timeseries to estimate delta

(g) alpha: significance level for the domain identification algorithm

(h) k: number of nearest neighbors to each grid cell i. The nearest neighbors are computed using the Haversine distance (https://en.wikipedia.org/wiki/Haversine_formula)

(i) tau_max: it defines the range of lags used in the network inference (i.e., for each pair of domains signals A and B, the code will test the statistical significance of correlations in the lag range R \in [-tau_max,tau_max])

(l) q: FDR parameter to test the significance of the lag-correlations (e.g., q = 0.05 implies that (on average) only 5% of the links identified is expected to be a false positive).

# (iv) OUTPUTS 

Outputs are saved in a "outputs" folder.
The folder "outputs" will contain 3 subfolders:

(a) seed_identification

Here you will find 

* local_homogeneity_field.npy

values of the local homogeneity of each grid cell i

* seed_positions.npy

The identified seeds: if a grid cell i, is not a seed it will be equal to 0. If a grid cell i is a seed, it will be equal to 1.

(b) domain_identification

Here you will find

* domain_ids.npy

a list with the domain ids

* domain_maps.npy

a numpy array with all the "domains map". E.g., if we identify N domains in a climate field with a x b grid points, the array defined in domain_maps will have the following format: np.shape(array) = [N,a,b].
Each entry array[i] will have the show the array with id defined in domain_ids[i]. 
Points belonging to a domain are equal to 1. 
In the folder "Notebooks" we show an example on how to work with such file.

(c) network_inference

Here you will find

* network_list.npy

It contains the list of edges identified.
Each edge have the following format:

[domain a, domain b, tau_min, tau_max, tau*, r*, weight]

where:

(i) r* is the max significant correlation

(ii) tau* is the lag correspondent to r*

(iii) [tau_min,tau_max] defines the range of lags associated with significant correlations
in the interval [r*-bartlett std(r*),r*+bartlett std(r*)]

(iv) if the range [tau_min, tau_max] include zero: domain a <---> domain b

if tau_min > 0:                               domain a  ---> domain b

if tau_max < 0:                               domain a <---  domain b

(vii) weight is the link weight: covariance at tau*

* strength_list.npy

list with domain strengths. Each entry in the list have the following format:
[domain id, strength]

* strength_map.npy

A map where each domain is defined by its strength.
In the folder "Notebooks" we show an example on how to plot the strength map.

# (v) PUBLICATIONS 

I. Fountalis, C. Dovrolis, A. Bracco, B. Dilkina, and S. Keilholz.δ-MAPS from spatio-temporal data to a weighted and lagged network between functional domain.Appl. Netw. Sci., 3:21, 2018.

Bracco, A.,Falasca,F., Nenes, A., Fountalis, I., and Dovrolis, C. (2018), Advancing climate science with Knowledge-Discovery through Data mining,NPJ Climate and Atmosph.Science,4, doi:10.1038/s41612-017-0006-4

Falasca,F., Bracco, A., Nenes, A., and Fountalis, I. (2019).  Dimensionality reduction and network inference for climate data usingδ-MAPS: Application to the CESM LargeEnsemble sea surface temperature. Journal of Advances in Modeling Earth Systems,11, 1479–1515. https://doi.org/10.1029/2019MS001654

Falasca,F., Crétat, J., Braconnot, and Bracco, A. Spatiotemporal complexity and time-dependent networks in sea surface temperature from mid- to late Holocene. Eur. Phys. J.Plus 135:392 (2020). https://doi.org/10.1140/epjp/s13360-020-00403-x
(free version here: https://www.researchgate.net/publication/341107349_Spatiotemporal_complexity_and_time-dependent_networks_in_sea_surface_temperature_from_mid-_to_late_Holocene)

