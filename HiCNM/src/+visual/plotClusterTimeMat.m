classdef plotClusterTimeMat < handle
    % ------------------------------------- %
    % --- plot cluster assignment      ---- %
    % ----@created 2013-09-24 EK       ---- %
    % --- @revised 2014-10-04 EK       ---- %
    % --- @revised 2021-09-09 DN       ---- %
    % ------------------------------------- %
    properties (Hidden)
        Name
        Data
        crossing
        Plot
    end
    methods
        function obj = plotClusterTimeMat(Data,identifier)           % Constructor
            obj.Name = 'CTM';
            obj.Name = [obj.Name,'_',identifier];
            obj.Data = Data;
            obj.Plot = 'Matrix';
        end
        
        function delete(obj)                     % Destructor
        end
        
        function Name = getName(obj)
            Name = obj.Name;
        end
        function setName(obj, Name)
            obj.Name = Name;
        end
        
        function run(obj,fig_handle)
            plotData(obj.Data, obj.Plot);
        end
    end
end


function plotData(Data,type)

% ------------------------------------- %
% --- plot cluster Time. matrix    ---- %
% ----@created 2020-09-09 Nan ---------- %
% ------------------------------------- %

%% Options/Parameters
% options:
ColorMap = getColormapForCTM();
ColorMap = vertcat([1 1 1],ColorMap);

%% Plot
switch type
    case 'Matrix'
        plotTimeMatrix(Data, ColorMap);
end

end