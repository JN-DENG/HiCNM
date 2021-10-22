function DTMC(CNM)

% ------------------------------------- %
% --- determine Transition Matrix ----- %
% ----@revised: 2021-09-03 ND check --- %
% ------------------------------------- %

%% Computes transition probability matrix and sorts clusters
[CNM.PM, CNM.c1_Labels, CNM.c1_Centroids] = determineClusterTransitionMat(CNM.ClusteringResults.c0_Labels, ...
    CNM.ClusteringResults.c0_Centroids,CNM.M,'basic', []);

%% Network, remove P_{ii}
for i=1:length(CNM.PM)
    if CNM.PM(i,i)~=1 
        CNM.PM(i,i)=0;
        CNM.PM(:,i)=CNM.PM(:,i)/sum(CNM.PM(:,i));
    end
end   

% for i=1:length(CNM.PM)
%     if CNM.PM(i,i)==1 % no back trans and make all 0
%         CNM.PM(:,i)=0;
%     else
%         CNM.PM(i,i)=0;
%         CNM.PM(:,i)=CNM.PM(:,i)/sum(CNM.PM(:,i));
%     end
% end  

%% remove swithching point of different trans
    slide_change = find( diff(CNM.c1_Labels)~=0)+1;
    label_trans=CNM.c1_Labels(slide_change);
    % Sequence of visiting clusters
    label_trans_new = [CNM.c1_Labels(1);label_trans];
    
    Num_swith=length(CNM.Data.t_switch);
    if Num_swith~=0
        N_fromswithing = zeros(size(CNM.PM));       
        N_all = zeros(size(CNM.PM));
        
        for i_switch=1:Num_swith
            trans_switch_point = CNM.Data.t_switch(i_switch)/CNM.Data.dt;
            trans_switch_point_previous = trans_switch_point-1;
            label_switch_point = CNM.c1_Labels(int32(trans_switch_point));
            label_switch_point_previous=CNM.c1_Labels(int32(trans_switch_point_previous));      
            ind=find(label_trans_new(2:end)==label_switch_point & label_trans_new(1:end-1)==label_switch_point_previous);
            N_all(label_switch_point,label_switch_point_previous)=length(ind);
            if isempty(ind) % No trans
            else % mulitiple trans 
                N_fromswithing(label_switch_point,label_switch_point_previous)=...
                    N_fromswithing(label_switch_point,label_switch_point_previous)+1; %remove this part  
            end 
        end
        P_correct=N_fromswithing./N_all;
        P_correct(isnan(P_correct))=0;
        
        % remove swithing probability
        CNM.PM=(1-P_correct).*CNM.PM;
        
        % normalize P
        for i_switch=1:Num_swith
            trans_switch_point_previous = trans_switch_point-1;
            label_switch_point_previous = CNM.c1_Labels(int32(trans_switch_point_previous));      
            if sum(CNM.PM(:,label_switch_point_previous))==0
                CNM.PM(:,label_switch_point_previous)=0;
            else
                CNM.PM(:,label_switch_point_previous)=CNM.PM(:,label_switch_point_previous)/...
                    sum(CNM.PM(:,label_switch_point_previous));        
            end
        end
    end
    
    
%% Network, Add P_{ii}=1 for permenant state
% for i=1:length(CNM.PM)
%     if sum(CNM.PM(:,i))==0
%         CNM.PM(i,i)=1;
%     end
% end 
    
%% Computes transition Time matrix and sorts clusters
[CNM.TimeM] = determineClusterTransitionTime(CNM.c1_Labels,CNM.Data.dt,CNM.Data.t_switch);
        
% Computes how sparse the transition matrix is
CNM.sparsity = determineSparsityOfCTM(CNM.PM);

% Cluster probability matrix qk = #snapshots in cluster k / #total snapshots
CNM.q = determineClusterProbVec(CNM.c1_Labels);
end