function [PMl] = determineDynamicsOfCTM(P, P_powers)
% ------------------------------------- %
% --- determine dynamics of CTM     --- %
% --- @created: 2013-09-27 EK --------- %
% --- @revised: 2013-10-11 EK --------- %
% --- @revised: 2014-02-11 EK rm type - %
% --- Piterl(:,:,l) = P^l;          --- %
% ------------------------------------- %

%% Dynamics of CTM
PMl = zeros(size(P,1),size(P,2),P_powers);
for l = 1:P_powers
    PMl(:,:,l) = P^l;
end

