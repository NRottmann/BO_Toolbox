classdef Branin
% 2-dmensional function: the sum of a quartic polynomial and a sinusoid of
% the first variable
% The range is x: [-5, 10] and y:[0, 15]
% Optima is around 0.3979 at [3.1416, 2.2750], [9.4248, 2.4750] and
% [−3.1416, 12.2750]
% reported in Michael James Sasena. Flexibility and Efficiency Enhancements for Constrained Global
% Design Optimization with Kriging Approximations. PhD thesis, University of Michigan,
% 2002.
%
% https://www.sfu.ca/~ssurjano/branin.html
%
% Date: 09.08.2019
% Author: Michael Werner   
    properties
        num_vars = 2;
        vars = [optimizableVariable('x',[-5,10]), optimizableVariable('y',[0,15])];
        minimize = true;
    end
    
    properties(Access=private)
       a = 1;
       b = 5.1/(4*pi^2);
       c = 5/pi;
       r = 6;
       s = 10;
       t = 1/(8*pi);
    end
    
    methods
        function obj = Branin()
           % 
        end
        
        function z = call(obj,x, y)
            %CALL calculates the value of the branin function at x and y
            z = obj.a*(y - obj.b*x^2 + obj.c*x - obj.r)^2 + ...
                obj.s*(1 - obj.t)*cos(x) + obj.s;
        end
    end
end