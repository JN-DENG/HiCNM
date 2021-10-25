%% Example: 4 transitions of fluidic Pinball at Re=80
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
%% Characteristic frequency
figure(1)
Fs = 10;  % Sampling Freequency 10Hz 
wdw = 1024;
[FFp,frq] = pwelch(detrend(a1(:,2),0),wdw,fix(0.95*wdw),wdw,Fs);
semilogy(frq,abs(FFp),'k','LineWidth',2);
grid on
xlim([0 1]);
xticks(0:0.2:1);
grid on
set(gca,'FontSize',20);
xlabel('$f$','Interpret','latex','FontSize',24);
ylabel('psd','Interpret','latex','FontSize',24);

Fc = frq(  FFp==max(FFp) );
disp(['Characteristic frequency: ',num2str(Fc)])

%% Filtering the data
close all

% Butterworth filter design 
n = 5;            % Filter order
Wn = Fc/Fs*(1/5); % Cutoff frequency
[bb,aa] = butter(n,Wn,'low'); % lowpass filter with cutoff frequency

Y1 = filtfilt(bb,aa,a1);
Y2 = filtfilt(bb,aa,a2);
Y3 = filtfilt(bb,aa,a3);
Y4 = filtfilt(bb,aa,a4);
Y=[Y1;Y2;Y3;Y4];

% Filtered POD amplitudes
a_f=[Y1;Y2;Y3;Y4];

% Input Data
dt = 0.1;
x = a_f;
t=0.1:dt:size(x,1)/10;
%% Prepare Input Data & options 
DatatoCNM.dt = dt;
DatatoCNM.t  = t;
DatatoCNM.ts = x;
DatatoCNM.t_switch=t([15001, 30001, 40001]);  % cut points for different transtions

%  Cluster analysis parameters
params_user.nClusters      	     = 20;         % Cluster number
params_user.nRepetitions         = 30;         % Iteration time
params_user.optimalClustering    = 'sparsity';

%  Transition matrix parameters
params_user.ClusterOrdering      ='transitions';

% Additional settings
params_user.save    = 1;
params_user.verbose = 0;
params_user.plot    = 1;

%% CROM
CNMobj = CNM(DatatoCNM,params_user);
CNMobj.run

save CNMobj_Base20.mat CNMobj
copyfile('output','L1')

