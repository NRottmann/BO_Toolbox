classdef HoelderTable
% Date: 09.08.2019
% Author: Michael Werner   
    properties
        num_vars = 2;
        vars = [optimizableVariable('x',[-10,10]), optimizableVariable('y',[-10,10])];
        minimize = true;
    end
    
    methods
        function obj = HoelderTable()
            %
        end
        
        function z = call(obj,x, y)
            %CALL calculates the value of the HoelderTable function at x and y
            z = - abs(sin(x) * cos(y) * exp(abs(1-sqrt(x^2*y^2)/pi)));
        end
    end
end