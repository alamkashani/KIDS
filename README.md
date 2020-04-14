# KIDS
A Fuzzy Density-Grid-Based Method for Clustering K-dimensional Data


we propose a novel fuzzy density-grid-based clustering method for clustering k-dimensional data based on the concept of diffusion. KIDS, an acronym for K-dimensional Ink Drop Spread, relies on discovering densely-connected pieces of data in k-dimensional grids, giving it the ability to simultaneously exploit the advantages of fuzziness, as well as both density-based and grid-based clustering. In the proposed method, each input data record is first mapped to a grid cell. The data points are then spread in k dimensions, just like ink drop patterns. Eventually, the impacts of all data grid cells are condensed and subjected to a threshold in order to compute the final clusters. This approach makes high-speed data clustering feasible without lowering the quality. The experimental results show that the proposed method has superior quality and efficiency for both low and high dimensions. In addition, the method is not only robust to noise, but is also capable of finding clusters of arbitrary shapes.

#Files:

kIDS_p2_Agg.m         :The KIDS method implemented on Aggregation dataset

kIDS_p2_Chainlink.m   :The KIDS method implemented on Chainlink dataset

Aggregation.txt       :The Aggregation dataset

chainlink.txt         :The Chainlink dataset
