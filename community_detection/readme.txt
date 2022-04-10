Finding spatiotemporal patterns with community detection. Here we use infomap.

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
################### Note:
#######################################################

(i) community detection algorithms will assign every grid cell to a community.
This is true even if some time series are actually noise.
(ii) communities do not have to be spatially contiguous.

d-MAPS take care of issues (i) and (ii).

If (i) and (ii) are not an issue for your spectific project, then infomap represents
a powerful was to reduce the dimensinality of your dataset.

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
