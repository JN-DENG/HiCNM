%% Example: 4 transitions of fluidic Pinball at Re=80
% Local dynamics in L2
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

c_label=CNMobj.c1_Labels;

clear CNMobj

for Index_SubCluster=1:max(c_label)
    a_Sub=a(c_label==Index_SubCluster,:);
    load(['loop_Base20_C',num2str(Index_SubCluster),'.mat'])

    %% PCA
    nClusters=size(CNMobj.c1_Centroids,1);
    c1_Centers=CNMobj.c1_Centroids;
    c1_Centers_mean = sum(CNMobj.c1_Centroids,1)./nClusters;
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
    a_r = a_Sub*V_r;
    a_r_f=CNMobj.Data.ts*V_r;
    %% FIGURE 
    close all
    Idx=CNMobj.c1_Labels;

    cmap=distinguishable_colors(nClusters);
    figure
    set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.8])
    set(gca,'position',[0.15,0.15,0.8,0.8] ); 
    %
    xlim([floor(min(a_r(:,2))) ceil(max(a_r(:,2)))])
    ylim([floor(min(a_r(:,1))) ceil(max(a_r(:,1)))])
    hold on
    % snapshots
    for i=1:nClusters 
        scatter(a_r(Idx==i,2),a_r(Idx==i,1),'.','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
        'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    end
    % Pos info
    posAxes = get(gca, 'Position');
    posX = posAxes(1);
    posY = posAxes(2);
    width = posAxes(3);
    height = posAxes(4);
    limX = get(gca, 'Xlim');
    limY = get(gca, 'Ylim');
    minX = limX(1);
    maxX = limX(2);
    minY = limY(1);
    maxY = limY(2);

    %arrow
    P=CNMobj.PM;
    %
    for i=1:10
        P(i,i)=0;
    end
    %
    [row,col]=find(P~=0);

    for i=1:length(row)
        p1=col(i);
        p2=row(i);
        if sum(sum([p1,p2]==[row,col],2)==2)==1
            % arc
            Center=(([Vec_r(p1,2),Vec_r(p1,1)]+...
                     [Vec_r(p2,2),Vec_r(p2,1)]))./2;
            R=norm( ([Vec_r(p2,2),Vec_r(p2,1)]-...
                     [Vec_r(p1,2),Vec_r(p1,1)]) );
            norm_vec=([Vec_r(p2,1),-Vec_r(p2,2)]-...
                [Vec_r(p1,1),-Vec_r(p1,2)])./R;
            xc=Center(1);
            yc=Center(2); 
            Center_rond=[xc + 0.5*sqrt(3)*R*norm_vec(1), ...
                yc + 0.5*sqrt(3)*R*norm_vec(2)];
            vec= [Vec_r(p1,2),Vec_r(p1,1)]-Center_rond;
            th1 = acos(dot(vec,[1,0])/(norm(vec)*norm([1,0])));
            if vec(2)<0
                th1=-th1;
            end
            th2 = th1-pi/3;
            th = linspace(th1,th2,101);
            x=Center_rond(1) + R*cos(th);
            y=Center_rond(2) + R*sin(th);
            plot(x,y,'color',[48 48 48]/255,...
                'LineWidth',6*P(p2,p1)); 
            xData = [Vec_r(p1,2),x(51)];
            yData = [Vec_r(p1,1),y(51)];
            x0 = xData(1 : 2);
            y0 = yData(1 : 2);
            x0(2)=x0(1)+1.05*(x0(2)-x0(1));
            y0(2)=y0(1)+1.05*(y0(2)-y0(1));
            xNew = posX + (x0 - minX) / (maxX - minX) * width;
            yNew = posY + (y0 - minY) / (maxY - minY) * height;
            annotation('arrow', xNew, yNew,'color',[48 48 48]/255,...
                'LineStyle','none',...
                'HeadLength', 15, 'HeadWidth', 18, 'HeadStyle', 'vback1')
        else
            xData = [Vec_r(p1,2),Vec_r(p2,2)];
            yData = [Vec_r(p1,1),Vec_r(p2,1)];
            x0 = xData(1 : 2);
            y0 = yData(1 : 2);
            x0(2)=x0(1)+0.5*(x0(2)-x0(1));
            y0(2)=y0(1)+0.5*(y0(2)-y0(1));
            xNew = posX + (x0 - minX) / (maxX - minX) * width;
            yNew = posY + (y0 - minY) / (maxY - minY) * height;
            plot(xData,yData,'color',[48 48 48]/255,...
                'LineWidth',6*P(p2,p1)); 
            annotation('arrow', xNew, yNew, 'color',[48 48 48]/255,...
                'LineStyle','none',...
                'HeadLength', 15, 'HeadWidth', 18, 'HeadStyle', 'vback1')
        end
    end

    %centroids
    for k = 1:nClusters  
        plot(Vec_r(k,2),Vec_r(k,1),'s','MarkerEdgeColor','k','MarkerFaceColor',cmap(k,:), 'MarkerSize',15)
    end
    for iCluster = 1:nClusters  
       text(Vec_r(iCluster,2)+0.06, ...
            Vec_r(iCluster,1)+0.04,num2str(iCluster),'FontSize',18);
    end

    set(gca,'FontSize',24);
    xlabel('$\gamma_2$','interpreter','latex','FontSize',32)
    ylabel('$\gamma_1$','interpreter','latex','FontSize',32)
    %%
    name = sprintf('%s%1d%s','./L2_C',Index_SubCluster,'/2Dmapping.png');
    print(1,name,'-dpng','-r300')
    % name = sprintf('%s','./2Dmapping_B20c1.eps');
    % print(1,name,'-depsc','-r300')
end
