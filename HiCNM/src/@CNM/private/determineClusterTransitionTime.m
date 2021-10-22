function [Time] = determineClusterTransitionTime(c0_Labels,dt,t_switch)
% ---------------------------------- %
% --- Cluster transition Time matrix ---- %
% ----@created: 2020-08-20 ND        %
% ----@revised: 2021-09-03 ND        %
% ---------------------------------- %

unique_labels = unique(c0_Labels);
nCluster      = length(unique_labels);

%% Initialization
Time  = zeros(nCluster,nCluster);          % cluster transition matrix

%% Cluster probability vector of c0_Labels
slide_change = find( diff(c0_Labels)~=0)+1;
label_trans=c0_Labels(slide_change);
%%first trans
slide_change_new=[1;slide_change];
label_trans_new = [c0_Labels(1);label_trans];

for jCluster = 1:nCluster
    for iCluster = 1:nCluster
        ind=find(label_trans_new(2:end)==iCluster & label_trans_new(1:end-1)==jCluster);
        if isempty(ind)
            Time(iCluster,jCluster)=0;
        else
            %t_ij = 0.5*(ti+tj) Transient time from j cluster to i cluster
%             if iCluster==label_trans(end)
%                 Time(iCluster,jCluster)=dt*0.5*mean(slide_change(ind(1:end-1)+2)-slide_change(ind(1:end-1)));
%             else
%                 Time(iCluster,jCluster)=dt*0.5*mean(slide_change(ind+2)-slide_change(ind));
%             end
            %t_ij = tj Residing time in the j cluster
            Time(iCluster,jCluster)=dt*mean(slide_change_new(ind+1)-slide_change_new(ind));
        end
    end
end

%% remove swithching point of different trans
Num_swith=length(t_switch);

    if Num_swith==0
    else 
        Time_fromswithing = zeros(size(Time));
        N_all = zeros(size(Time));
        N_fromswithing = zeros(size(Time));

        for i_switch=1:Num_swith
            trans_switch_point = t_switch(i_switch)/dt;
            trans_switch_point_previous = trans_switch_point-1;
            
            label_switch_point=c0_Labels(int32(trans_switch_point));
            label_switch_point_previous=c0_Labels(int32(trans_switch_point_previous)); 
         
            ind=find(label_trans_new(2:end)==label_switch_point & ...
                     label_trans_new(1:end-1)==label_switch_point_previous);
            N_all(label_switch_point,label_switch_point_previous)=length(ind);
    
            ind_finalstate = find(slide_change_new==int32(trans_switch_point)) -1;
            
            NUM_finalstate=slide_change_new(ind_finalstate+1)-slide_change_new(ind_finalstate);
            
            if isempty(ind)
            else
                Time_fromswithing(label_switch_point,label_switch_point_previous)=...
                    Time_fromswithing(label_switch_point,label_switch_point_previous)+...
                    dt*NUM_finalstate;
                N_fromswithing(label_switch_point,label_switch_point_previous)=...
                    N_fromswithing(label_switch_point,label_switch_point_previous)+1; %remove this part             
            end
        end
        
        Time_correct = zeros(size(Time));

        for jCluster = 1:nCluster
            for iCluster = 1:nCluster
                if N_fromswithing(iCluster,jCluster)==0
                    Time_correct(iCluster,jCluster)=Time(iCluster,jCluster); 
                else
                    if N_all(iCluster,jCluster)-N_fromswithing(iCluster,jCluster) ==0
                        Time_correct(iCluster,jCluster)=0;
                    else
                        Time_correct(iCluster,jCluster)=...
                            (Time(iCluster,jCluster)*N_all(iCluster,jCluster)...
                            -Time_fromswithing(iCluster,jCluster))/...
                            (N_all(iCluster,jCluster)-N_fromswithing(iCluster,jCluster));
                    end 
                end
            end
        end
        Time=Time_correct;
    end

end