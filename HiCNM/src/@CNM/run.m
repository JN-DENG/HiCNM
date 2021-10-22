function run(CNM, varargin)
% run CNM pipeline
% PURPOSE : CNM execution program which executes the main functions and some analysis stuff
% INPUT   : CNM : object of type CNM
% NOTE    : If mode save == 1, writes results into mat files

if isempty(varargin) == 1
    % nix
end

if exist(utils.Parameters.instance.parameters.path2save) ~= 7
    mkdir(utils.Parameters.instance.parameters.path2save);
end

ClusterAnalysis(CNM);

Analysis(CNM);

if utils.Parameters.instance.parameters.plot == 1
    plotResults(CNM);
end
display(CNM);

end