classdef Ackley
% The Ackley benchmark function
%
% Date: 15.08.2019
% Author: Michael Werner   

    properties(Access=public)
        num_vars = 2;
        vars = [optimizableVariable('x',[-5,5]), optimizableVariable('y',[-5,5])];
        minimize = true;
    end
    
    methods
        function obj = Ackley()
            %
        end
        
        function s = call(obj, x, y)
            %CALL calculates the value of the Ackley function          
            s = - 20 * exp(-0.2*sqrt(0.5*(x^2+y^2)))...
                - exp(0.5*(cos(2*pi*x) + cos(2*pi*y)))...
                + exp(1) + 20;
            
            if obj.minimize
                s = -s;
            end
        end
    end
end