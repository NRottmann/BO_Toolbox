classdef LinearFeature
    %LinearFeature Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=public)
        num_param
    end
    
    properties(Access=private)
       vars
       D
       d
       T
    end
    
    methods(Access=public)
        function obj = LinearFeature(D, d, inputVars, varargin)
            %NEURALNET Construct an instance of this class
            %   Detailed explanation goes here
            
            % Default values
            if d ~= 1
               error("LinearFeature can only provide a one dimensional feature space"); 
            end
            defaultargs = {}; 
            params = setargs(defaultargs, varargin);
            obj.D = D;
            obj.d = d;
            obj.T = ones(1, obj.D) * 1/sqrt(obj.D);
            obj.vars = inputVars;
            obj.num_param = 0;
        end
               
        function f = getfeature(obj, x, param)
           if len(param) ~= obj.num_param
               error(strcat("The given number of parameters does not match the requested size of ", num2str(obj.num_param), "."));
           end          
           f = obj.T*x;
        end
               
        function bounds = getbounds(obj)
            %GETBOUNDS returns the bounds of the feature space           
            % activation function 1/(1+exp(-a)) returns values from 0 to 1
            bounds = zeros(1, 2);
            for l=1:obj.D
                bounds(1,1) = bounds(1,1) + obj.vars(l).Range(1) * obj.T(1,l);
                bounds(1,2) = bounds(1,2) + obj.vars(l).Range(2) * obj.T(1,l);
            end
        end
    end
end