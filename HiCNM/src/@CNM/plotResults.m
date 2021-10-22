function plotResults(CNM)

%% Cluster index
Data.Data = CNM.c1_Labels;
Data.Time = CNM.Data.t;
Data.T    = max(CNM.Data.t);
Data.T_switch = CNM.Data.t_switch-CNM.Data.dt;
plotThisClass = visual.plotLabels(Data); 
CNM.plotRunner(plotThisClass);

%% Cluster transition matrix
plotThisClass = visual.plotClusterProbMat(CNM.PM, 'P'); % crossing possible
CNM.plotRunner(plotThisClass);

%% Cluster transition time matrix
plotThisClass = visual.plotClusterTimeMat(CNM.TimeM, 'Time'); % crossing possible
CNM.plotRunner(plotThisClass);

%% Cluster centroids distance matrix
% plotThisClass = visual.plotClusterMatrix(CNM.DM); % crossing possible
% CNM.plotRunner(plotThisClass);

%% Voronoi diagram
% clear Data
% Data.data       = CNM.ts_r;
% Data.labels     = CNM.c1_Labels;
% Data.centroids  = CNM.c1_Centroids_r;
% plotThisClass   = visual.plotVoronoiDiagram(Data);
% CNM.plotRunner(plotThisClass);

%% Evolution of CTM
% for i = 1:size(CNM.PMl,3)
%     plotThisClass = visual.plotClusterProbMat(CNM.PMl(:,:,i), ['Pl_',sprintf('%02g',i)],[]); % crossing possible
%     CNM.plotRunner(plotThisClass);
% end

%% Evolution of probability distribution
% Data.data = CNM.pl;
% Data.time = [0:length(CNM.pl{1})];
% plotThisClass = visual.plotEvolutionOfProbability(Data);
% for iCluster = 1:length(CNM.pl)
%      plotThisClass.setIndex(iCluster);
%      CNM.plotRunner(plotThisClass);
% end
end