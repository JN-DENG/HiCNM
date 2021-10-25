%% Example: 4 transitions of fluidic Pinball at Re=80
%% RUN in just one click

% Initialize and clean
clear;
clc; close all;
%% Clustering analysis

% HiCNM in L1: Clustering analysis with filtered data
L1_Base20_filtered

% HiCNM in L2: Clustering analysis with filtered data
L2_SubCluster_loop

%% Visualization
% Directed graph and centroids: L1
L1_Directed_Graph
L1_Plot_centroids

% Directed graph and centroids: L2
L2_Directed_Graph
L2_Plot_centroids
