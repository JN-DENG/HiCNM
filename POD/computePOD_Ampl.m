%% Initialize and clean
clear;
clc; close all;

% load node and element information
addpath('../Dataset/');
Grid2=load('Mesh/Grid2.dat','-ascii');
elem=load('Mesh/elem.dat','-ascii');
coef=load('Mesh/find.dat','-ascii');

NDE = length(Grid2);  % number of grids
NEL = length(elem)/4; % number of T6 triangular elements
Nxy = length(elem);   % number of sub-elements
dt =1; % time stepping with a number of snapshot 
%% grid elements coordinates
coord1=Grid2(elem(:,1),:);
coord2=Grid2(elem(:,2),:);
coord3=Grid2(elem(:,3),:);
c=(coord1 + coord2 + coord3)./3; % Centroid of elements

% compute the area of each triangles
A = zeros(Nxy,1);
for j=1:Nxy
    mat=[coord2(j,1)-coord1(j,1),coord3(j,1)-coord1(j,1);coord2(j,2)-coord1(j,2),coord3(j,2)-coord1(j,2)];
    A(j)=0.5*abs(det(mat));
end

% Base  flow (U,V)_0; here is the symmetric steady solution
ss_SYM = load('SS_SYM', '-ascii');
us = ss_SYM(:,1);
vs = ss_SYM(:,2);

% POD Modes
load POD_UALL_symmetrized.mat
load POD_VALL_symmetrized.mat

%% --- Projection to the leading 400 POD modes SYM---

disp('---------------- Start Calculating for SYM Trans. -------------------')
tic

load UALL_tran_Re80_sym.mat
load VALL_tran_Re80_sym.mat

% interpolate velocity fields at the center of the triangles
uPODm = (U_POD(elem(:,1),:)+U_POD(elem(:,2),:)+U_POD(elem(:,3),:))/3; % velocity in center
vPODm = (V_POD(elem(:,1),:)+V_POD(elem(:,2),:)+V_POD(elem(:,3),:))/3; % velocity in center

Udif = UALL(:,1:dt:end) - us;
Vdif = VALL(:,1:dt:end) - vs;
    
% define velocity at the center of the triangles
Udifm = (Udif(elem(:,1),:)+Udif(elem(:,2),:)+Udif(elem(:,3),:))/3;
Vdifm = (Vdif(elem(:,1),:)+Vdif(elem(:,2),:)+Vdif(elem(:,3),:))/3;
    
a = ( (A.*uPODm)'*Udifm + (A.*vPODm)'*Vdifm )';
save POD_Ampl_Re80_sym.dat a -ascii;

% plot leading 20 POD mode amplitudes
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i=1:20
  subplot(5,4,i)
  plot(a(:,i))
  xlabel('time')
  ylabel(['a_{',num2str(i),'}'])
  title(['Mode Amplitude of POD No ',num2str(i)])
end
saveas(1,'POD_Ampl_SYM.png')

disp('Symmetry Transient data -- COMPLETE!')
toc


%% --- Projection to the leading 400 POD modes UP---
close all;
clear UALL VALL
disp('---------------- Start Calculating for ASYM Trans. from the UPwards deflected steady solution -------------------')
tic

load UALL_tran_ASYM80_UP.mat
load VALL_tran_ASYM80_UP.mat

% interpolate velocity fields at the center of the triangles
uPODm = (U_POD(elem(:,1),:)+U_POD(elem(:,2),:)+U_POD(elem(:,3),:))/3; % velocity in center
vPODm = (V_POD(elem(:,1),:)+V_POD(elem(:,2),:)+V_POD(elem(:,3),:))/3; % velocity in center

Udif = UALL(:,1:dt:end) - us;
Vdif = VALL(:,1:dt:end) - vs;
    
% define velocity at the center of the triangles
Udifm = (Udif(elem(:,1),:)+Udif(elem(:,2),:)+Udif(elem(:,3),:))/3;
Vdifm = (Vdif(elem(:,1),:)+Vdif(elem(:,2),:)+Vdif(elem(:,3),:))/3;
    
a = ( (A.*uPODm)'*Udifm + (A.*vPODm)'*Vdifm )';
save POD_Ampl_Re80_up.dat a -ascii;

% plot leading 20 POD mode amplitudes
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i=1:20
  subplot(5,4,i)
  plot(a(:,i))
  xlabel('time')
  ylabel(['a_{',num2str(i),'}'])
  title(['Mode Amplitude of POD No ',num2str(i)])
end
saveas(1,'POD_Ampl_up.png')

disp('UP Transient data -- COMPLETE!')
toc


%% --- Projection to the leading 400 POD modes DOWN---
close all;
clear UALL VALL
disp('---------------- Start Calculating for ASYM Trans. from the DOWNwards deflected steady solution -------------------')
tic

load UALL_tran_ASYM80_DOWN.mat
load VALL_tran_ASYM80_DOWN.mat

% interpolate velocity fields at the center of the triangles
uPODm = (U_POD(elem(:,1),:)+U_POD(elem(:,2),:)+U_POD(elem(:,3),:))/3; % velocity in center
vPODm = (V_POD(elem(:,1),:)+V_POD(elem(:,2),:)+V_POD(elem(:,3),:))/3; % velocity in center

Udif = UALL(:,1:dt:end) - us;
Vdif = VALL(:,1:dt:end) - vs;
    
% define velocity at the center of the triangles
Udifm = (Udif(elem(:,1),:)+Udif(elem(:,2),:)+Udif(elem(:,3),:))/3;
Vdifm = (Vdif(elem(:,1),:)+Vdif(elem(:,2),:)+Vdif(elem(:,3),:))/3;
    
a = ( (A.*uPODm)'*Udifm + (A.*vPODm)'*Vdifm )';
save POD_Ampl_Re80_down.dat a -ascii;

% plot leading 20 POD mode amplitudes
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i=1:20
  subplot(5,4,i)
  plot(a(:,i))
  xlabel('time')
  ylabel(['a_{',num2str(i),'}'])
  title(['Mode Amplitude of POD No ',num2str(i)])
end
saveas(1,'POD_Ampl_down.png')

disp('DOWN Transient data -- COMPLETE!')
toc

%% --- Projection to the leading 400 POD modes SYM mirror ---
close all;
clear UALL VALL
disp('---------------- Start Calculating for SYM mirror Trans. -------------------')
tic

load UALL_tran_Re80_sym.mat
load VALL_tran_Re80_sym.mat

% mirror
Nsnap=size(UALL,2); % Number of Snapshots
UmALL=zeros(size(UALL));
VmALL=zeros(size(UALL));

tic
%  LOOP ON THE FRAMES   
   for i=1:Nsnap 
       Umir = zeros(NDE,1);
       Vmir = zeros(NDE,1);
%  Velocity EXTRACTION from the flow file
       Uin = UALL(:,i);
	     Vin = VALL(:,i);
        for j=1:8633
            Umir(j)=   Uin(coef(j,1))*coef(j,4)+Uin(coef(j,2))*coef(j,5)+Uin(coef(j,3))*coef(j,6);
            Vmir(j)= -(Vin(coef(j,1))*coef(j,4)+Vin(coef(j,2))*coef(j,5)+Vin(coef(j,3))*coef(j,6));
        end
        UmALL(:,i)=Umir;
        VmALL(:,i)=Vmir;
   end
% steady solution
usmir = zeros(NDE,1);
vsmir = zeros(NDE,1);
for j=1:8633
    usmir(j)=   us(coef(j,1))*coef(j,4)+us(coef(j,2))*coef(j,5)+us(coef(j,3))*coef(j,6);
    vsmir(j)= -(vs(coef(j,1))*coef(j,4)+vs(coef(j,2))*coef(j,5)+vs(coef(j,3))*coef(j,6));
end
disp('-------------------- mirror COMPLETE ---------------------')  
toc

%
clear UALL VALL
UALL=UmALL;
VALL=VmALL;

% interpolate velocity fields at the center of the triangles
uPODm = (U_POD(elem(:,1),:)+U_POD(elem(:,2),:)+U_POD(elem(:,3),:))/3; % velocity in center
vPODm = (V_POD(elem(:,1),:)+V_POD(elem(:,2),:)+V_POD(elem(:,3),:))/3; % velocity in center

Udif = UALL(:,1:dt:end) - us;
Vdif = VALL(:,1:dt:end) - vs;
    
% define velocity at the center of the triangles
Udifm = (Udif(elem(:,1),:)+Udif(elem(:,2),:)+Udif(elem(:,3),:))/3;
Vdifm = (Vdif(elem(:,1),:)+Vdif(elem(:,2),:)+Vdif(elem(:,3),:))/3;
    
a = ( (A.*uPODm)'*Udifm + (A.*vPODm)'*Vdifm )';
save POD_Ampl_Re80_mirror_sym.dat a -ascii;

% plot leading 20 POD mode amplitudes
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i=1:20
  subplot(5,4,i)
  plot(a(:,i))
  xlabel('time')
  ylabel(['a_{',num2str(i),'}'])
  title(['Mode Amplitude of POD No ',num2str(i)])
end
saveas(1,'POD_Ampl_mirror_SYM.png')

disp('SYM Mirror Transient data -- COMPLETE!')
toc
