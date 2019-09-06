classdef PCAFeature
    %NEURALNET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=public)
        num_param = 0;
    end
    
    properties(Access=private)
        d
    end
    
    methods(Access = private)
       
    end
    
    methods(Access=public)
        function obj = PCAFeature(D, d, inputVars, varargin)
            %PCAFeature Construct an instance of this class
            %   Detailed explanation goes here
            obj.d = d;
        end
               
        function f = getfeature(obj, x, param)
            x_in = x;
            % augment data when the number of observations is smaller than
            % the number of variables
            x = x';
            while size(x, 1) <= size(x,2)
                x = [x; x];
            end
            
            % fit pca
            [coeffs] = pca(x);
            % transform
            f = x_in' * coeffs;
            % extract d first components
            f = f(:, 1:obj.d);
            % scale into range of [0, 1]
            f = (f - min(f(:))) / (max(f(:)) - min(f(:)));
            f = f';
        end
               
        function bounds = getbounds(obj)
            %GETBOUNDS returns the bounds of the feature space           
            bounds = [zeros(obj.d, 1), ones(obj.d, 1)];
        end
    end
end