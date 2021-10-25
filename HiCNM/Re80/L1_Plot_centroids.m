%% Example: 4 transitions of fluidic Pinball at Re=80
% Initialize and clean
clear;
clc; close all;

addpath('../src/');

% Prepare Data & options 
DatatoCNM.dt = 1;
DatatoCNM.t  = 1;
DatatoCNM.ts = 1;

%  Cluster analysis parameters
params_user.nClusters      	     = 20;
params_user.nRepetitions         = 30;
params_user.optimalClustering    = 'sparsity';
params_user.rDim                 = 3;
params_user.rVec                 = [1,2,3];

%  Transition matrix parameters
params_user.ClusterOrdering      ='transitions';

% Additional settings
params_user.save    = 1;
params_user.verbose = 0;
params_user.plot    = 1;

%% CNM load with struct
CNMobj = CNM(DatatoCNM,params_user);

load CNMobj_Base20.mat

figure
plot(CNMobj.c1_Labels,'.')
grid on

%% centroid in real space
a1=load('../../POD/POD_Ampl_Re80_sym.dat');
a2=load('../../POD/POD_Ampl_Re80_mirror_sym.dat');
a3=load('../../POD/POD_Ampl_Re80_down.dat');
a4=load('../../POD/POD_Ampl_Re80_up.dat');
%
a1=a1(1:15000,:); % 15000 snapshots for symmetry trajectory
a2=a2(1:15000,:); % 15000 snapshots for symmetry mirrored trajectory
a3=a3(1:10000,:); % 10000 snapshots for asymmetry down trajectory
a4=a4(1:10000,:); % 10000 snapshots for asymmetry up trajectory

a=[a1;a2;a3;a4];

% Centroids by averaging the original data
c1_Centroids=zeros(20,400);
for i=1:20
    c1_Centroids(i,:)=sum( a(CNMobj.c1_Labels==i,:),1)/sum(CNMobj.c1_Labels==i);
end

%% steady solution
ss_SYM = load('../../Dataset/SS_SYM', '-ascii');
us = ss_SYM(:,1);
vs = ss_SYM(:,2);

%% CENTROIDS
load ../../POD/POD_UALL_symmetrized.mat
load ../../POD/POD_VALL_symmetrized.mat

% Centroids from Clustering algorithm (filtered data)
% UALL_Centroids=us+U_POD(:,1:400)*(CNMobj.c1_Centroids)';
% VALL_Centroids=vs+V_POD(:,1:400)*(CNMobj.c1_Centroids)';

% Centroids by averaging the original data
UALL_Centroids=us+U_POD(:,1:400)*(c1_Centroids)';
VALL_Centroids=vs+V_POD(:,1:400)*(c1_Centroids)';
Nt=size(UALL_Centroids,2);
% load node and element information
load ../../Dataset/Mesh/Grid2.dat
load ../../Dataset/Mesh/elem.dat

%% Vorticity calculation
disp('Velocity field loaded,')
disp('Calculating the vorticity field ...')
NDE = 8633;
NEL = 4225;

X=Grid2(:,1);
Y=Grid2(:,2);
CONN=elem; % CONNECTIVITY MATRIX

% calculate vorticity field
VORTALL = zeros (NDE,Nt);

for i=1:Nt
    U=UALL_Centroids(:,i);
    V=VALL_Centroids(:,i);

    DX=0 ;
    DY=0 ;

    VORT = zeros (NDE,1);
    IND = zeros (NDE,1);
    
    for KC=1:4*NEL
        K1= CONN(KC,1);
        K2= CONN(KC,2);
        K3= CONN(KC,3);
        K4= CONN(KC,1);
        
        X1 = X(K1);
        X2 = X(K2);
        X3 = X(K3);
        X4 = X(K4);
        DX = max(DX, max([X1,X2,X3,X4])-min([X1;X2;X3;X4]));
        Y1 = Y(K1);
        Y2 = Y(K2);
        Y3 = Y(K3);
        Y4 = Y(K4);
        DY = max(DY, max([Y1;Y2;Y3;Y4])-min([Y1;Y2;Y3;Y4]));
        U1 = U(K1);
        V1 = V(K1);
        U2 = U(K2);
        V2 = V(K2);
        U3 = U(K3);
        V3 = V(K3);
        U4 = U(K4);
        V4 = V(K4);

        AREA = ( (X3-X1)*(Y4-Y2) + (X2-X4)*(Y3-Y1) ) ;
        CIRC = ( (U1+U2)*(X2-X1) + (V1+V2)*(Y2-Y1)...
            + (U2+U3)*(X3-X2) + (V2+V3)*(Y3-Y2)...
            + (U3+U4)*(X4-X3) + (V3+V4)*(Y4-Y3)...
            + (U4+U1)*(X1-X4) + (V4+V1)*(Y1-Y4) ) / AREA ;

        VORT(K1) = VORT(K1) + CIRC ;
        IND(K1) = IND(K1) + 1;
        VORT(K2) = VORT(K2) + CIRC ;
        IND(K2) = IND(K2) + 1;
        VORT(K3) = VORT(K3) + CIRC ;       
        IND(K3) = IND(K3) + 1;
        if K4 ~= K1 
            VORT(K4) = VORT(K4) + CIRC ;
            IND(K4) = IND(K4) + 1;
        end

    end
    VORTALL(:,i)=VORT./IND;
end
disp('Vorticity field  --  Completed!')


%% Visualization
% close all;
J=redblueTecplot;
disp('Visualization...')
Nstep = 1 ; % time stepping
figure(2)
for i=1:Nstep:20
    subplot(4,5,i)
    h = trisurf (CONN, X, Y, VORTALL(:,i), 'facecolor','interp') ;
    set(gca,'DataAspectRatio',[1 1 1]);
    view(2) ;
    axis tight; 
    %colorbar
    caxis([ -1.5 1.5 ]) ;
    shading interp
    colormap(J)
    %
    title(['Centroid', num2str( i )])
end
% saveas(2,'Centroids.png')

%%
rondphi=0:pi/100:2*pi+pi/100;
R=0.5;
rondx=R*cos(rondphi);
rondy=R*sin(rondphi);

for i=1:20
    close all
    figure
    hold all
    set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.35])
    set (gca,'position',[0.05,0.05,0.9,0.9] );  
    h = trisurf (CONN, X, Y, VORTALL(:,i), 'facecolor','interp');
    set(gca,'DataAspectRatio',[1 1 1]);
    view(2) ;
    axis tight; 
    %colorbar
    caxis([ -1.5 1.5]) ;
    shading interp
    colormap(J)
    xlim([-4,20])
    ylim([-4,4])
    set (gca,'xtick',[]); 
    set (gca,'ytick',[]); 
    % rond 
    fill(rondx,0.75+rondy,'k')
    fill(rondx,-0.75+rondy,'k')
    fill(-0.75*sqrt(3)+rondx,0+rondy,'k')
    % borders of regions
    for iRegions = 1:4
        plot3([-4,-4], [-4,4],[1,1],  '-k', 'LineWidth',1) % horizontal up
        plot3([20,20], [-4,4],[1,1],  '-k', 'LineWidth',1) % horizontal below
        plot3([-4,20], [-4,-4],[1,1], '-k','LineWidth',1.5) % vertical right
        plot3([-4,20], [4,4],[1,1],   '-k', 'LineWidth',1.5) % vertical left
    end
    hold off
    % save
    pic_name=['./L1/Centoirds_',num2str(i),'.png'];
    saveas(1,pic_name)
end
