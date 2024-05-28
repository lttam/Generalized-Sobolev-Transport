% This is the code for the paper
%       "Generalized Sobolev Transport for Probability Measures on a Graph"
%       ICML2024
%-------------------------------------------------------------------

% Third-party toolbox
% + Sobolev transport
% + Sinkhorn_log

% Using the third-party toolbox to create the graphs G_Log and G_Sqrt
% amazon_OT_1000_RandLLE.mat (G_Log)
% amazon_OT_1000_RandSLE.mat (G_Sqrt)
% For examples:
% -- For building random graph (G_Log / G_Sqrt) from support data points
% (Sobolev transport third-party toolbox)
% + clusteringDataset_buildRandomGraph_Log: build random connected graph G_Log (M
% nodes, and M log(M) edges)
% + clusteringDataset_buildRandomGraph_Sqrt: build random connected graph
% G_Sqrt (M nodes, and M^(3/2) edges)

% Compute_GST_EXP1_vec_V3: compute GST with phi1 = @(X) exp(X) - X - 1;
% Compute_GST_EXP2_vec_V3: compute GST with phi2 = @(X) exp(X.^2) - 1;
% OrliczWasserstein: compute Orlicz-Wasserstein by Guha et al., 2023
% (Algorithm 1)

% Samples
% + test_gst1_all.m : for GST with phi1
% + test_gst2_all.m : for GST with phi2

% + test_OrliczWasserstein_Phi1.m : for OW with phi1
% + test_OrliczWasserstein_Phi2.m : for OW with phi2

% Data:
% amazon_ID10000.mat: the 10K random pairs

% Time consumption for GST, OW with phi1, phi2
% + time_GenGraphSobolev_EXP1_V3_1T.m : GST-phi1
% + time_GenGraphSobolev_EXP1_V3_1T.m : GST-phi2
% + time_OW_EXP1_GraphMetric.m : OW-phi1
% + time_OW_EXP1_GraphMetric.m : OW-phi2

% For the limit case: phi0
% GST is equal to Sobolev transport (ST)
% OW is equal to OT with graph metric as its ground cost

% -- Note:
% + The code uses Graph and Network Algorithms from MATLAB. (e.g., Dijkstra
% algorithm for shortest path from a source point to a destination set of
% points.
% + The code uses Optimization Toolbox from MATLAB. (e.g., fmincon to solve
% the univariate optimization problem for GST).



% ================================================================
% DROPBOX shared link for 'amazon_1000_RandLLE_Graph.mat' (~510MB)
% https://www.dropbox.com/scl/fi/osuv7on4bw4p8pijmemf4/amazon_1000_RandLLE_Graph.mat?rlkey=796udajwf0177gwsbhdwvjohz&dl=0




