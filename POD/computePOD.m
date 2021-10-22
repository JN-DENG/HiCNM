% Initialize and clean
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
%% (U,V) Velocity field of 3 transitions

% % velocity field for transition starting with the SYMMETRIC steady solution
load UALL_tran_Re80_sym.mat
load VALL_tran_Re80_sym.mat
UALL1 = UALL(:,1:10:15000);
VALL1 = VALL(:,1:10:15000);

% % starting with the UPwards deflected steady solution
load UALL_tran_ASYM80_UP.mat
load VALL_tran_ASYM80_UP.mat
UALL2 = UALL(:,1:10:10000);
VALL2 = VALL(:,1:10:10000);

% % starting with the DOWNwards deflected steady solution
load UALL_tran_ASYM80_DOWN.mat
load VALL_tran_ASYM80_DOWN.mat
UALL3 = UALL(:,1:10:10000);
VALL3 = VALL(:,1:10:10000);
%% Base  flow (U,V)_0; here is the symmetric steady solution
ss_SYM = load('SS_SYM', '-ascii');
us = ss_SYM(:,1);
vs = ss_SYM(:,2);

%% mirror reflect the data for transition starting with the SYMMETRIC steady solution 
Nx=size(UALL1,1);
Nm=size(UALL1,2);
UmALL=zeros(size(UALL1));
VmALL=zeros(size(UALL1));
tic
%  LOOP ON THE SNAPSHOTS   
   for i=1:Nm 
       Umir = zeros(Nx,1);
       Vmir = zeros(Nx,1);
%  Velocity EXTRACTION from the flow file
       Uin = UALL1(:,i);
	   Vin = VALL1(:,i);
        for j=1:8633
            Umir(j)=   Uin(coef(j,1))*coef(j,4)+Uin(coef(j,2))*coef(j,5)+Uin(coef(j,3))*coef(j,6);
            Vmir(j)= -(Vin(coef(j,1))*coef(j,4)+Vin(coef(j,2))*coef(j,5)+Vin(coef(j,3))*coef(j,6));
        end
        UmALL(:,i)=Umir;
        VmALL(:,i)=Vmir;
   end
% steady solution
usmir = zeros(Nx,1);
vsmir = zeros(Nx,1);
for j=1:8633
    usmir(j)=   us(coef(j,1))*coef(j,4)+us(coef(j,2))*coef(j,5)+us(coef(j,3))*coef(j,6);
    vsmir(j)= -(vs(coef(j,1))*coef(j,4)+vs(coef(j,2))*coef(j,5)+vs(coef(j,3))*coef(j,6));
end
   
disp('-------------------- Mirror data COMPLETE ---------------------')  
toc

%% fluctuating velocity field
Uprime = [UALL1- us,UmALL-usmir,UALL2- us,UALL3- us];
Vprime = [VALL1- vs,VmALL-vsmir,VALL2- vs,VALL3- vs];
Nsnap=size(Uprime,2); % Number of Snapshots
clear UALL1 UALL2 UALL3 UmALL VALL1 VALL2 VALL3 VmALL 

% grid elements coordinates
coord1=Grid2(elem(:,1),:);
coord2=Grid2(elem(:,2),:);
coord3=Grid2(elem(:,3),:);
c=(coord1 + coord2 + coord3)./3; % Centroid of elements

% interpolate velocity fields at the center of the triangles
um = (Uprime(elem(:,1),:)+Uprime(elem(:,2),:)+Uprime(elem(:,3),:))/3; % velocity in center
vm = (Vprime(elem(:,1),:)+Vprime(elem(:,2),:)+Vprime(elem(:,3),:))/3;
% compute the area of each triangles
A = zeros(Nxy,1);
for j=1:Nxy
    mat=[coord2(j,1)-coord1(j,1),coord3(j,1)-coord1(j,1);coord2(j,2)-coord1(j,2),coord3(j,2)-coord1(j,2)];
    A(j)=0.5*abs(det(mat));
end

%% POD to compress the data set
% Correlation matrix
C = ( (A.*um)'*um + (A.*vm)'*vm )/Nsnap;
%%compute POD after subtracting mean (i.e., do PCA)
[PSI,S] = eig(C);
% [V,D] = eig(A), A*V = V*D. V eigenvector, D eigenvalue
lambda = diag(real(S));

% Proportion of energy of POD modes
figure
semilogy(lambda./sum(lambda),'LineWidth',2)
hold on
semilogy(lambda./sum(lambda),'.r','MarkerSize',12)
xlabel('Mode number')
ylabel('\lambda')
xlim([0,100])

% Proportion of energy of the leading POD modes
Energy=zeros(1,Nsnap);
for i=1:Nsnap
    Energy(i)=sum(lambda(1:i))/sum(lambda);
end
%
figure
plot(1:Nsnap,Energy,'LineWidth',2);
hold on
plot(1:Nsnap,Energy,'.r','MarkerSize',12);
xlabel('Mode number')
ylabel('Energy percent')
xlim([0,100])
ylim([.99,1])

% POD modes u_i = 1/sqrt(M*lambd_i) * sum_{m=1}^M(V_i^m U_prime^m)
U_POD = real( (Uprime*PSI)./(sqrt(Nsnap*lambda))' );
V_POD = real( (Vprime*PSI)./(sqrt(Nsnap*lambda))' );
% Mode amplitudes a_i = sqrt(M*lambd_i) * V_i^m
a = real(PSI)*diag(sqrt(lambda*Nsnap));

disp('Proportion of energy of the leading 400 POD modes: ')
sprintf('%.16f',Energy(400))

% Vorticity computation
VORT_POD = Comp_Vorticity(U_POD(:,1:20),V_POD(:,1:20),Grid2,elem);

%% Data sauvage
% save('CorrelationMatrix.mat','C')

% save leading 400 POD modes
U_POD=U_POD(:,1:400);
save('POD_UALL_symmetrized.mat','U_POD')
V_POD=V_POD(:,1:400);
save('POD_VALL_symmetrized.mat','V_POD')

%% Plot the POD modes
X=Grid2(:,1);
Y=Grid2(:,2);
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
vortmin = -1;  % only plot what is in -0.1 to 0.1 range
vortmax = 1;
for i=1:20
   subplot(5,4,i)
   VORT=VORT_POD(:,i);
   VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
   VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin
   h = trisurf (elem, X, Y, VORT, 'facecolor','interp','EdgeColor','none');
   set(gca,'DataAspectRatio',[1 1 1]);
   view(2) ;
   axis tight; 
   colorbar;
   caxis([-0.3 0.3])
   shading interp 
   colormap ('jet');
   title(['Vorticity field of POD No ',num2str(i)])   
end
saveas(2,'POD_leading 20 modes.png')
