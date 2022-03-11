# py-dMAPS

---

Contacts

- Fabrizio Falasca (fabrifalasca@gmail.com) and Ilias Fountalis (Foudalisi@hotmail.com)

Info on limitations/bugs and general possible issues about the code are in the folder "debugging_etc". 
The folder is managed by Ljuba Novi (ljubanovi@yahoo.it)

---

δ-MAPS is a method that identifies the semi-autonomous components (i.e., patterns) of a spatiotemporal climate field and studies their weighted and potentially lagged interactions. The components of the system are modeled as spatially contiguous, possibly overlapping, homogeneous domains. To the best of our knowledge, δ-MAPS is the first method that can identify spatially contiguous and possibly overlapping clusters (i.e., domains). At a second step, δ-MAPS infers a directed and weighted network between the identified domains. Edge direction captures the temporal ordering of events while the weight associated to each edge captures the magnitude of interaction between domains.

The framework can be used also as a preprocessing step for running causal algorithms on climate time series. Causal graph discovery often benefits from dimensionality reduction. Domains signals (i.e., time series) identified by δ-MAPS can then be fed to causal discovery algorithms (e.g., https://jakobrunge.github.io/tigramite/) to identify causal linkages.

---

# Contents

(i) Required Python Packages

(ii) Netcdf files 

(iii) Preprocessing

(iv) how to run py-dMAPS

(v) Outputs

(vi) Tutorials

(vii) Updates

(viii) Publications

# (i) Required Python packages

(a) numpy

(b) scipy

(c) scikit-learn

(d) netCDF4



# (ii) Netcdf files 

The code accepts netcdf files https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_introduction.html .
Each time series x_i(t) in the dataset is assumed to be already preprocessed (i.e., each x_i(t) should be stationary). In our case we work often with monthly anomalies: the preprocessing involves (a) trend removal and (b) removal of the seasonal cycle.

To see the structure of a netcdf file "file.nc", open a terminal and type: 

ncdump -h file.nc

To visualize netcdf files please consider the software NCVIEW by David Pierce (freely available at http://meteora.ucsd.edu/~pierce/ncview_home_page.html). Another useful software is Panoply by NASA (freely available at https://www.giss.nasa.gov/tools/panoply/).

# (iii) PREPROCESSING 

To preprocess a netcdf dataset we suggest to use the CDO package from MPI (https://code.mpimet.mpg.de/projects/cdo/).

Below we show how to preprocess a spatiotemporal dataset saved as monthly averages.
Consider a spatiotemporal climate field such as the COBEv2 reanalysis (https://psl.noaa.gov/data/gridded/data.cobe2.html) with temporal and spatial resolution of 1 month and 1 by 1 degree (i.e., 360x180 grid points) respectively. The temporal range goes from January 1850 to December 2018.

Possible preprocessing steps using CDO:

* Remap the dataset to the resolution you want (i.e, 2 by 2 degree)

cdo -L -remapbil,r180x90 input.nc output.nc

* Select a period (e.g., from January 1950 to December 2010)

cdo selyear,1950/2010 input.nc output.nc

* Remove the seasonal cycle

cdo -L -ymonsub input.nc -ymonmean input.nc output.nc 

* Remove a linear trend

cdo detrend input.nc output.nc

* Focus only on the latitudinal range from 60 degree South to 60 degree North

cdo sellonlatbox,0,360,-60,60 input.nc output.nc

# (iv) HOW TO RUN py-dMAPS 

* python3 run_delta_maps.py -i configs/sample_config.json

Inputs in the configs/sample_config.json

(a) path_to_data: path to the spatiotemporal dataset. Each time series x_i(t) in the dataset is assumed to be already preprocessed (i.e., each x_i(t) should be stationary)

(b) output_directory: name of the output directory that will be created

(c) variable name: name of the variable of interest in the netcdf file (i.e., "sst" in our case)

(d) latitude: name of the variable containing latitudes (i.e., lat in our case, but other dataset could have "latitude")

(e) longitude: name of the variable containing longitude (i.e., lon in our case, but other dataset could have "longitude")

(f) delta_rand_samples: random sample of pairs of timeseries to estimate delta

(g) alpha: significance level for the domain identification algorithm

(h) k: number of nearest neighbors to each grid cell i. The nearest neighbors are computed using the Haversine distance (https://en.wikipedia.org/wiki/Haversine_formula)

(i) tau_max: it defines the range of lags used in the network inference (i.e., for each pair of domains signals A and B, the code will test the statistical significance of correlations in the lag range R \in [-tau_max,tau_max])

(l) q: FDR parameter to test the significance of the lag-correlations (e.g., q = 0.05 implies that (on average) only 5% of the links identified is expected to be a false positive).

<!--
## ----------------- OBSERVATION/NEW FEATURE -----------------

The starting point for the domain identification step is the selection of cores. We say that a cell "i" is a core if its local homogeneity (i.e., average pairwise correlation between time series in the K-neighborhood of i, i included) is (a) a local maximum in its K-neighborhood and (b) it is greater than the threshold delta. Main reason for this choice is to be able to differentiate between what is noise and significant variability of the system (i.e., dMaps will not include all time series in a domain).

A different option is to start the domain identification using all grid cells as core. If this option is chosen, then ALL points will be assigned to a domain. This can be sometimes preferred as (a) it allows not to worry about sensitivity of cores on the K parameters and (b) in certain cases it is preferred to assign all grid cell to a certain partition.

If this is your choice for domain identification do the following:

(a) open SeedIdentification/seed_identification.py 

(b) go to line 102-103 and comment the following line:

if(seed_homogeneity > max(cell_homogeneities) and seed_homogeneity > delta):

    seed_positions[index_i,index_j] = 1;
    
(c) simply replace with: 

seed_positions[index_i,index_j] = 1;

In this way all grid cells will be assigned to a domain. Note: the algorithm will be WAY slower, so do this with a relatively low resolution (e.g., 180x30 points would take x minutes) (Fab: insert x when you have it).

-->

# (v) OUTPUTS 

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

# (vi) Tutorials

Here you will find tutorials on how to analyze the data.

(a) Test 1. 

Output of domains identification and network inference for COBEv2 reanalysis. Global sea surface temperature (SST) from January 1980 to December 2015. Temporal resolution: 1 month. Number of grid points: 180x60 (2 by 2 degree with no High lats). 

It took 4843.15 seconds to complete the run (domain identification + network inference).

Note: 
- For the COBEv2 dataset, the "true" number of grid points is almost 1/2 of 180x60 as we are focusing on SST and the land is masked. 
- Domains can overlap: when we plot the domains and strength maps in the tutorial it is not possible to see the overlapping. A way to show it is to plot also the domains borders (we will do it soon).
- Main point: d-MAPS tries to discover large (constrained to your data) scale modes of variability. Find the "minimum" resolution that take into account of the processes you are interested in and use that. E.g., if your interested is in discovering modes at global scale in SST a resolution of 2 by 2 degrees is probably enough for a good qualitative answer (going to 0.5 by 0.5 degrees would make things exponentially slower and you are probably going to find similar large-scale answers). 
- As the files.nc are > 25MB, we couldn't upload them in this github repository.

- How you choose the parameters is domain' specific. A general starting point is K = 4 (or 8), alpha = 0.01 and FDR parameter q = 0.05.

# (viii) PUBLICATIONS 

I. Fountalis, C. Dovrolis, A. Bracco, B. Dilkina, and S. Keilholz.δ-MAPS from spatio-temporal data to a weighted and lagged network between functional domain.Appl. Netw. Sci., 3:21, 2018.

Bracco, A.,Falasca,F., Nenes, A., Fountalis, I., and Dovrolis, C. (2018), Advancing climate science with Knowledge-Discovery through Data mining,NPJ Climate and Atmosph.Science,4, doi:10.1038/s41612-017-0006-4

Falasca,F., Bracco, A., Nenes, A., and Fountalis, I. (2019).  Dimensionality reduction and network inference for climate data usingδ-MAPS: Application to the CESM LargeEnsemble sea surface temperature. Journal of Advances in Modeling Earth Systems,11, 1479–1515. https://doi.org/10.1029/2019MS001654

Falasca,F., Crétat, J., Braconnot, and Bracco, A. Spatiotemporal complexity and time-dependent networks in sea surface temperature from mid- to late Holocene. Eur. Phys. J.Plus 135:392 (2020). https://doi.org/10.1140/epjp/s13360-020-00403-x
(free version here: https://www.researchgate.net/publication/341107349_Spatiotemporal_complexity_and_time-dependent_networks_in_sea_surface_temperature_from_mid-_to_late_Holocene)

Novi, L., Bracco, A. & Falasca, F. Uncovering marine connectivity through sea surface temperature. Sci Rep 11, 8839 (2021). https://doi.org/10.1038/s41598-021-87711-z

