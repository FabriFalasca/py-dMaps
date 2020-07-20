Notebook file with test for the net inference

The input of the code are 
(a) a file.txt containing the cumulative anomalies (signals) of domains previously identified by delta-MAPS
(b) the ids of these domains

The outputs are

(a) a network defined as a list.
Each sublist inside the list is formatted as follows:
[domain A, domain B, tau min, tau max, tau*, r*, weight]

where 
- domains A and B are 2 domains (duh...)
- r* is their maximum significant correlation
- tau* is the lag for which we have r*
- tau min and tau max define the range of lags for which we have (i) significant correlations and (ii) each correlations is NOT LARGER THAN r* + its bartlett standard deviation and NOT SMALLER THAN r* - its bartlett standard deviation
- weight is the covariance at that fucking tau*

(b) strength_list
    a list with the the strength of each fucking domain.
