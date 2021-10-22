function ClusterAnalysis(CNM)
% Call the clustering algorithm

ClusterRun  = cell(1,utils.Parameters.instance.parameters.nRepetitions);
kmeans_eval = zeros(utils.Parameters.instance.parameters.nRepetitions,1);
for iRun = 1:utils.Parameters.instance.parameters.nRepetitions
    % k-means
    [CNM.ClusteringResults.c0_Labels, CNM.ClusteringResults.c0_Centroids, ...
        CNM.ClusteringResults.sumD, CNM.ClusteringResults.D] = kmeans(CNM.Data.ts,...
        utils.Parameters.instance.parameters.nClusters, ...
        'MaxIter', utils.Parameters.instance.parameters.nIterations, ...
        'Distance', utils.Parameters.instance.parameters.distmetric, ...
        'Start', 'plus');
     
    DTMC(CNM);  
    
    %%
    ClusterRun{iRun}.P                  = CNM.PM;
    ClusterRun{iRun}.Time               = CNM.TimeM;
    ClusterRun{iRun}.c1_Labels          = CNM.c1_Labels;
    ClusterRun{iRun}.c1_Centroids       = CNM.c1_Centroids;
    ClusterRun{iRun}.q                  = CNM.q;
    ClusterRun{iRun}.sparsity           = CNM.sparsity;
    ClusterRun{iRun}.ClusteringResults  = CNM.ClusteringResults;
    
    if strcmp(utils.Parameters.instance.parameters.optimalClustering,'sparsity')
        kmeans_eval(iRun) = CNM.sparsity;
    else
        kmeans_eval(iRun)  = sum(CNM.ClusteringResults.sumD);
    end
end
[~,best_result] = max(kmeans_eval);

CNM.PM                     = ClusterRun{best_result}.P;
CNM.TimeM                  = ClusterRun{best_result}.Time;
CNM.c1_Labels              = ClusterRun{best_result}.c1_Labels;
CNM.c1_Centroids           = ClusterRun{best_result}.c1_Centroids;
CNM.q                      = ClusterRun{best_result}.q;
CNM.sparsity               = ClusterRun{best_result}.sparsity;
CNM.ClusteringResults      = ClusterRun{best_result}.ClusteringResults;

%% label sorted
AppearOrder=unique(CNM.c1_Labels,'stable'); %appearing order 
mapClusterIndexAO(AppearOrder)=1:utils.Parameters.instance.parameters.nClusters;
CNM.PM = CNM.PM(AppearOrder,AppearOrder);
CNM.TimeM = CNM.TimeM(AppearOrder,AppearOrder);
CNM.c1_Labels=mapClusterIndexAO(CNM.c1_Labels);
CNM.c1_Centroids= CNM.c1_Centroids(AppearOrder,:);
CNM.q = CNM.q (AppearOrder);
CNM.ClusteringResults.c0_Labels=mapClusterIndexAO(CNM.ClusteringResults.c0_Labels);
CNM.ClusteringResults.c0_Centroids = CNM.ClusteringResults.c0_Centroids(AppearOrder,:);
CNM.ClusteringResults.sumD = CNM.ClusteringResults.sumD(AppearOrder);
CNM.ClusteringResults.D = CNM.ClusteringResults.D(:,AppearOrder);
end