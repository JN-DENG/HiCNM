%% Example: 4 transitions of fluidic Pinball at Re=80
% Feature space 3D for checking the baseline flow
% Initialize and clean
clear;
clc; close all;

addpath('../src/');
path2figs = './output/'; mkdir(path2figs)

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
%% Prepare Data & options 
DatatoCNM.dt = dt;
DatatoCNM.t  = t;
DatatoCNM.ts = x;

%  Cluster analysis parameters
params_user.nClusters      	     = 30;
params_user.nRepetitions         = 10;
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
%% PCA
nClusters=size(CNMobj.c1_Centroids,1);
%%
% c1_Centers=CNMobj.c1_Centroids;
c1_Centers=zeros(30,400);
for i=1:20
    c1_Centers(i,:)=sum(a(CNMobj.c1_Labels==i,:),1)/sum(CNMobj.c1_Labels==i);
end
%%
c1_Centers_mean = sum(c1_Centers,1)./nClusters;
c1_Centers_fluct=zeros(size(c1_Centers));
for iCluster = 1:nClusters
    c1_Centers_fluct(iCluster,:) = c1_Centers(iCluster, :) - c1_Centers_mean;
end
B = c1_Centers_fluct' * c1_Centers_fluct;


[V,Lambda]  = eig(B);
lambda      = diag(Lambda);
[lambda,IX] = sort(lambda,'descend');
V           = V(:,IX);

% Step 4 - Representation in R^r: first r principal components
nDim=3;
lambda_r     = lambda(1:nDim);
V_r          = V(:,1:nDim);

% Step 5 - Projection of centroids onto principal components
c1_Centers_r = c1_Centers * V_r;

% Projection of ai
Vec_r = c1_Centers_r;
a_r = a*V_r;
a_r_f = CNMobj.Data.ts*V_r;
%% FIGURE 
close all
Idx=CNMobj.c1_Labels;

cmap=distinguishable_colors(nClusters);
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.8])
set (gca,'position',[0.15,0.15,0.8,0.8] ); 
%
xlim([floor(min(a_r(:,1))) ceil(max(a_r(:,1)))])
ylim([floor(min(a_r(:,2))) ceil(max(a_r(:,2)))])
zlim([floor(min(a_r(:,3))) ceil(max(a_r(:,3)))])
%
hold on
% snapshots
for i=1:nClusters 
    scatter3(a_r(Idx==i,1),a_r(Idx==i,2),a_r(Idx==i,3),'.','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
    'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
end

%centroids
for k = 1:nClusters  
    plot3(Vec_r(k,1),Vec_r(k,2),Vec_r(k,3),'s','MarkerEdgeColor','k','MarkerFaceColor',cmap(k,:), 'MarkerSize',15)
end
hold off
view(45,30)
set(gca,'FontSize',24);
xlabel('$\gamma_1$','interpreter','latex','FontSize',32)
ylabel('$\gamma_2$','interpreter','latex','FontSize',32)
zlabel('$\gamma_3$','interpreter','latex','FontSize',32)

%%
% name = sprintf('%s','./FeatureSpace_3D.png');
% print(1,name,'-dpng','-r300')
