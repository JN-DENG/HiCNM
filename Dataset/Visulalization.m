% Initialize and clean
clear;
clc; close all;

% load node and element information
addpath('./Mesh/');
load Grid2.dat
load elem.dat
%% 3 transitions available at Re=80
% % velocity field for transition starting with the asymmetric steady solution
% % DOWNwards deflected
% load UALL_tran_ASYM80_DOWN.mat
% load VALL_tran_ASYM80_DOWN.mat
% t_start=0;

% % UPwards deflected
% load UALL_tran_ASYM80_UP.mat
% load VALL_tran_ASYM80_UP.mat
% t_start=0;

% % velocity field for transition starting with the SYMMETRIC steady solution
load UALL_tran_Re80_sym.mat
load VALL_tran_Re80_sym.mat
t_start=0;

%% 3 steady solutions

% % DOWNwards deflected
% load SS_ASYM_DOWN

% % UPwards deflected
% load SS_ASYM_UP

% % SYMMETRIC steady solution
% load SS_SYM

%% Vorticity calculation
disp('Velocity field loaded,')
disp('Calculating the vorticity field ...')
Nt = length(UALL);      % time step
NDE = length(Grid2);  % number of grids
NEL = length(elem)/4; % number of T6 triangular elements

X = Grid2(:,1); Y = Grid2(:,2);

% Snapshot to check
index_snapshot = 15000; % snapshot number
U = UALL(:,index_snapshot);
V = VALL(:,index_snapshot);

% calculate the vorticity field
VORT = Comp_Vorticity(U,V,Grid2,elem);

%% Visulalization
close all;
J=redblueTecplot;
disp('Visulalization...')

% cylinder
rondphi=0:pi/100:2*pi+pi/100;
R=0.5;
rondx=R*cos(rondphi);
rondy=R*sin(rondphi);

% 
figure
hold on
set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.35])
set (gca,'position',[0.05,0.05,0.9,0.9] );  
h = trisurf (elem, X, Y, VORT, 'facecolor','interp');
set(gca,'DataAspectRatio',[1 1 1]);
view(2) ;
axis tight; 
% colorbar
caxis([ -1.5 1.5]);
shading interp
colormap(J)
xlim([-4,20])
ylim([-4,4])
set (gca,'xtick',[]); 
set (gca,'ytick',[]); 
% cylinders 
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

%% Video for all the snapshots
    
% % calculate the vorticity field for all the snapshots (take long time and large workspace)
% VORTALL = Comp_Vorticity(UALL,VALL,Grid2,elem);
% disp('Vorticity field  --  Completed!')
% save('VORTALL.mat','VORTALL')
% disp('Vorticity field saved')   

% Video making
close all
figure
% create the video writer with 1 fps
writerObj = VideoWriter('Video_Re80_sym.avi');
writerObj.FrameRate = 10; % set the seconds per image
% open the video writer
open(writerObj);   

Nstep = 1000 ; % time stepping for a number of snapshots
for i=1:Nstep:Nt
    clf
    %
    U=UALL(:,i);
    V=VALL(:,i);
    VORT=Comp_Vorticity(U,V,Grid2,elem);
    hold on
    h = trisurf (elem, X, Y, VORT, 'facecolor','interp','EdgeColor','none') ;
    set(gca,'DataAspectRatio',[1 1 1]);
    view(2) ;
    axis tight; 
    colorbar;
    caxis([ -1.5  1.5 ]) ;
    shading interp
    colormap(J)
    % cylinders 
    fill(rondx,0.75+rondy,'k')
    fill(rondx,-0.75+rondy,'k')
    fill(-0.75*sqrt(3)+rondx,0+rondy,'k')
    hold off
    %
    title(['Vorticity field for time = ', num2str( 0.1*(i-1) )])
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    drawnow
end
% close the writer object
close(writerObj);


