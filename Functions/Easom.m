classdef Easom
% Date: 09.08.2019
% Author: Michael Werner   
    properties
        num_vars = 2;
        vars = [optimizableVariable('x',[-100,100]), optimizableVariable('y',[-100,100])];
        minimize = true;
    end
    
    methods
        function obj = Easom()
            
        end
        
        function z = call(obj,x, y)
            %CALL calculates the value of the Easom function at x and y
            z = -cos(x)*cos(y)*exp(-((x-pi)^2+(y-pi)^1));
        end
    end
end