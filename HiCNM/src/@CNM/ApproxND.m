function ApproxND(CNM, rDim, iAlpha)
if nargin == 1
    rDim  = [2]; %utils.Parameters.instance.parameters.rDim
    iAlpha = [1 2]; %utils.Parameters.instance.parameters.rVec
end

[CNM.ts_r,CNM.c1_Centroids_r,CNM.pca_vec_r] = ...
    CNM.compLowOrderRepresentation(CNM.Data.ts,CNM.c1_Centroids,rDim, iAlpha);

end