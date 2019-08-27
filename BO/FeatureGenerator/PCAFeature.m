classdef PCAFeature
    %NEURALNET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=public)
        num_param % TODO
    end
    
    properties(Access=private)

    end
    
    methods(Access = private)
       
    end
    
    methods(Access=public)
        function obj = PCAFeature(D, d, inputVars, varargin)
            %PCAFeature Construct an instance of this class
            %   Detailed explanation goes here
           % TODO 
        end
               
        function f = getfeature(obj, x, param)
           % TODO
        end
               
        function bounds = getbounds(obj)
            %GETBOUNDS returns the bounds of the feature space           
            % activation function 1/(1+exp(-a)) returns values from 0 to 1
            bounds = 0 % TODO
        end
    end
end