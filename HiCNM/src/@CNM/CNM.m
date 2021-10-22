classdef CNM < handle
    % CNM Class
    % INPUT : {filename} : m-file script with parameters
    %          Dataobject: object with data
    properties (SetAccess = public, GetAccess = public)
        Data              % time-series data
        M                 % Snapshots number
        ClusteringResults % Kmeans result before sorting by order of appearance  
        c1_Labels         % Cluster index of time series after sorting
        c1_Centroids      % Centroids
        q                 % data distribution in clusters
        PM                % Cluster transition probability matrix
        TimeM             % Cluster transition time matrix
        DM                % Cluster distance matrix
        sparsity          % P matrix sparsity
        ts_r
        c1_Centroids_r
        pca_vec_r
        PMl
        pl
    end
    
    properties (Dependent = true, SetAccess = private)
    end
    
    methods (Static)
         plotRunner(plotClass);
    end
    
    methods (Hidden)
        function checkProperties(obj)
            keyboard
        end
    end
    
    methods (Access = private)
        setParameters(params_or_filename);
    end
    
    methods(Static)
        [ts_r,c1_Centroids_r,pca_vec_r] = compLowOrderRepresentation(ts,c1_Centroids,rDim,rVec);
    end
    methods
        run(obj);                              % Run CNM pipeline
        ClusterAnalysis(obj);
        DTMC(obj);
        ApproxND(obj);
        GeometricProperties(obj);
        Dynamics(obj);
        Analysis(obj);                          % Contains all calls for the analysis
        display(obj);
        plotResults(obj);
        
        function obj = CNM(varargin)           % Constructor
            if not(isempty(varargin))
                if length(varargin) == 2
                    setParameters(varargin{2});
                    obj.Data = varargin{1};
                    obj.M    = size(obj.Data.ts,1);
                elseif length(varargin) == 1
                    if ischar(varargin{1}) == 1
                        setParameters(varargin{1});
                    else
                        obj.Data = varargin{1};
                        obj.M    = size(obj.Data.ts,1);
                    end
                end
            end
            
        end
        
        function delete(obj)                     % Destructor
        end
        
        function clean(obj)
            if exist(utils.Parameters.instance.parameters.path2save) == 7 % is folder
                rmdir(utils.Parameters.instance.parameters.path2save,'s');
                mkdir(utils.Parameters.instance.parameters.path2save);
            else
                mkdir(utils.Parameters.instance.parameters.path2save);
            end
        end
        
        function saveCNM(obj,filename)
            if nargin == 2
%                 file_CNM = ['CNM_results_',filename];
%                 file_params = ['CNM_parameters_',filename];
                file_CNM = [filename,'_results'];
                file_params = [filename,'_parameters'];
            else
                file_CNM = ['CNM_results'];
                file_params = ['CNM_parameters'];
            end
            save(file_CNM, 'obj');
            params = utils.Parameters.instance.parameters;
            save(file_params,'params');
        end
        
        function [CNM,Parameters]= loadCNM(obj,name1, name2)
            if nargin < 2
                name1 = 'CNM_results';
                name2 = 'CNM_parameters';
            end
            CNMstruct = load(name1);
            CNM = getfield(CNMstruct, 'obj');
            
            Params_struct = load(name2);
            Parameters = getfield(Params_struct, 'params');
            
            p = utils.Parameters.instance();
            p.parameters = Parameters;
            %obj.
        end
        
    end
    
end


