function Dynamics(CNM)

% Powers of the transitionmatrix
P_powers = utils.Parameters.instance.parameters.powerLmat;
[CNM.PMl] = determineDynamicsOfCTM(CNM.PM, P_powers);

% Evolution of probability vector for different initial conditions
L_powers = utils.Parameters.instance.parameters.powerL;
nCluster = utils.Parameters.instance.parameters.nClusters;
CNM.pl  = cell(nCluster,1);
for iCluster = 1:nCluster
    CNM.pl{iCluster,1} = determineDynamicsOfSPV(iCluster, CNM.PM, L_powers);
end

% % Asymptotic probability vector
% CNM.pinf = CNM.P^(10^4)*[1/nCluster.*ones(nCluster,1)];

end