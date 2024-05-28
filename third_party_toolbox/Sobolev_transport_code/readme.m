% This is the code for the paper 
% "Sobolev transport: a scalable metric for probability measures with graph metrics"
% AISTATS 2022

% Data: e.g., 'twitter.mat' for the TWITTER dataset.

% Third-party toolbox
% + figtreeKCenterClustering: for the farthest point clustering (used for quantization supports into M clusters)
% + mexEMD: for optimal transport distance computation.

% -- For building random graph (G_Log / G_Sqrt) from support data points
% + clusteringDataset_buildRandomGraph_Log: build random connected graph G_Log (M
% nodes, and M log(M) edges)
% + clusteringDataset_buildRandomGraph_Sqrt: build random connected graph
% G_Sqrt (M nodes, and M^(3/2) edges)

% -- Compute distance matrices for Sobolev transport, optimal transport and
% tree-Wasserstein
% + compute_SobolevTransport: compute the distance matrix for Sobolev
% transport
% + compute_OT_GroundGraphMetric: compute the distance matrix for optimal
% transport with ground graph metric
% + compute_TW: compute the distance matrix for tree-Wasserstein with
% random tree metric extracted from the graph structure.

% -- Note:
% The code uses Graph and Network Algorithms from MATLAB. (e.g., Dijkstra
% algorithm for shortest path from a source point to a destination set of
% points.

