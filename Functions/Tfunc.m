classdef Tfunc
% 2-dmensional test function for Bayesian Optimization which has only one
% effective dimension
%
% Date: 10.09.2019
% Author: Michael Werner   

    properties(Access=public)
        num_vars = 2;
        vars = [optimizableVariable('x', [-10,10]),...
                optimizableVariable('y', [-10,10])];
        minimize = false;
    end
    
    properties(Access=private)
       T = [0.7071, 0.7071];
       a = [0.1 1 1];
    end
    
    methods
        function obj = Tfunc()

        end
        
        function z = call(obj, varargin)
            %CALL calculates the value of the TFunc function 
            input = cell2mat(varargin);
            v = [input(1); input(2)];
            
            val = obj.T * v; 

            z = -obj.a(1)*val^2 + obj.a(2)*val + obj.a(3);
        end
    end
end