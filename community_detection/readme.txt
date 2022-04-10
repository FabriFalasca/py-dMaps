Finding spatiotemporal patterns with community detection. We do so through the infomap method.
Community detection algorithms are an additional way to reduce the dimensionality of spatiotemporal data 
without using PCA or any other matrix decomposition algorithm. The are useful as they 
allow to focus on correlations rather than variance, they define patterns that are not orthogonal.

Therefore, even components with smaller variance are identified by the algorithm. This is useful, especially for nonlinear dynamics.

Additionally, if needed, results of d-MAPS can be compared to the one from community detection.

# Note

(i) community detection algorithms will assign every grid cell to a community. This is true even if some time series are actually noise.
(ii) communities do not have to be spatially contiguous.

d-MAPS take care of issues (i) and (ii).

If (i) and (ii) are not an issue for your specific project, then infomap is a powerful way to reduce the dimensionality of your dataset.

This readme is a work in progress.

#######################################################
################### How does it work:
#######################################################

(a) We start from a spatiotemporal dataset: time series embedded in a 2 dimensional
grid. The code creates a network between time series. 2 time series are connected
if their correlation is greater than a threshold K.

(b) K is then similar to delta in d-Maps and we use the same heuristics to find it.
In fact the only parameter here is alpha, significance level for the
delta (now K) heuristics.

(c) Once you have a network representation of your dataset (adjacency matrix)
we reduce it dimensionality via community detection. We use a powerful community
detection algorithm known as infomap to find communities.
Those communities are the patterns we are trying to find. And they are then
remapped back into the map.

#######################################################
################### Main packages needed
#######################################################

- networkx
https://github.com/networkx/networkx

- infomap
Give a look here: https://www.mapequation.org/infomap/


#######################################################
################### Some papers on infomap
#######################################################

- Paper by Tantet and Dijkstra using infomap for SST
https://esd.copernicus.org/articles/5/1/2014/

- Benchmarks for community detection algorithms showing very good
results with infomap

https://arxiv.org/abs/0908.1062

- Main infomap paper: https://www.mapequation.org/assets/publications/RosvallBergstromPNAS2008Full.pdf
Infomap website: https://www.mapequation.org/publications.html
