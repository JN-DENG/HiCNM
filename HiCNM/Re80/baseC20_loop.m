%% Example: 4 transitions of fluidic Pinball at Re=80
% Initialize and clean
clear;
clc; close all;

addpath('../src/');
path2figs = './output/'; 
[status, message, messageid] = rmdir(path2figs);
mkdir(path2figs)

% POD amplitudes loading
a1=load('../../POD/POD_Ampl_Re80_sym.dat');
a2=load('../../POD/POD_Ampl_Re80_mirror_sym.dat');
a3=load('../../POD/POD_Ampl_Re80_down.dat');
a4=load('../../POD/POD_Ampl_Re80_up.dat');

a1=a1(1:15000,:); % 15000 snapshots for symmetry trajectory
a2=a2(1:15000,:); % 15000 snapshots for symmetry mirrored trajectory
a3=a3(1:10000,:); % 10000 snapshots for asymmetry down trajectory
a4=a4(1:10000,:); % 10000 snapshots for asymmetry up trajectory

a=[a1;a2;a3;a4];

% Input Data
dt = 0.1;
x = a;
t=0.1:dt:size(x,1)/10;

%% Prepare Input Data & options 
DatatoCNM.dt = dt;
DatatoCNM.t  = t;
DatatoCNM.ts = x;
DatatoCNM.t_switch=t([15001, 30001, 40001]);  % cut points for different transtions

%  Cluster analysis parameters
params_user.nClusters      	     = 20;
params_user.nRepetitions         = 5;
params_user.optimalClustering    = 'sparsity';

%  Transition matrix parameters
params_user.ClusterOrdering      ='transitions';

% Additional settings
params_user.save    = 1;
params_user.verbose = 0;
params_user.plot    = 1;

%% CNM
CNMobj = CNM(DatatoCNM,params_user);
load CNMobj_Base20.mat

c_label = CNMobj.c1_Labels; % Evolution of cluster index
components = CNMobj.c1_Centroids; % Centroids

clear CNMobj

%% Loop for subclusters
dt = 0.1;
t0 = 0.1*(1:size(x,1));

for i=1:20
    disp(['Inlet of cluster ',num2str(i)])
    close all
    t_cluster_a=find( c_label==i);
    cut_point=find(diff(t_cluster_a)~=1)+1;
    
    cluster_a=a(t_cluster_a,:);
    x = cluster_a;
    t=0.1:dt:size(x,1)/10;
    figure(1)
    for j=1:20
        subplot(5,4,j)
        plot(t0,a(:,j))
        hold on
        plot(t0(t_cluster_a),a(t_cluster_a,j),'.k')
        plot(dt*t_cluster_a(cut_point),x(cut_point,j),'or')
        axis tight
    end
%     saveas(1,['output/Input_fluc_C',num2str(i),'.png'])
    
    %% Prepare Data & options 
    DatatoCNM.dt = dt;
    DatatoCNM.t  = t;
    DatatoCNM.ts = x;
    DatatoCNM.t_switch=t(cut_point);  % starting points of different transtion

    %  Cluster analysis parameters
    params_user.nClusters      	     = 10;
    params_user.nRepetitions         = 100;
    params_user.optimalClustering    = 'sparsity';

    %  Transition matrix parameters
    params_user.ClusterOrdering      ='transitions';

    % Additional settings
    params_user.save    = 1;
    params_user.verbose = 0;
    params_user.plot    = 1;

    %% CNM
    CNMobj = CNM(DatatoCNM,params_user);
    CNMobj.run
    
    filename=['loop_BaseC20_C',num2str(i),'.mat'];
    save(filename,'CNMobj')

    %%
    clustercenter=CNMobj.c1_Centroids;
    Idx=CNMobj.c1_Labels;
    figure(2);
    plot3(a(:,2),a(:,3),a(:,1),'.k')
    hold on;
    for j=1:params_user.nClusters  
    plot3(x(Idx==j,2),x(Idx==j,3),x(Idx==j,1),'.');
    end
    cVec = 'bgrcmykbgrcmyk';
    for k = 1:params_user.nClusters  
        plot3(clustercenter(k,2),clustercenter(k,3),clustercenter(k,1),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
    end
    xlim([min(x(:,2))-0.01,max(x(:,2)+0.01)])
    ylim([min(x(:,3))-0.01,max(x(:,3)+0.01)])
    zlim([min(x(:,1))-0.01,max(x(:,1)+0.01)])
%     saveas(2,['output/centroids_C',num2str(i),'.png'])
    copyfile('output',['L2_C_',num2str(i)])
end