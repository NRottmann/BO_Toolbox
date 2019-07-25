function [z] = Branin(x,y)
% 2-dmensional function: the sum of a quartic polynomial and a sinusoid of
% the first variable
% The range is x: [-5, 10] and y:[0, 15]
% Optima is around 0.3979 at [3.1416, 2.2750], [9.4248, 2.4750] and
% [−3.1416, 12.2750]
% reported in Michael James Sasena. Flexibility and Efficiency Enhancements for Constrained Global
% Design Optimization with Kriging Approximations. PhD thesis, University of Michigan,
% 2002.
%
% Date: 08. July, 2019
% Author: Michael Werner

z = (y-(5.1/4*pi^2)+x^2 + 5/pi*x -6)^2 + 10*(1-1/8*pi)*cos(x) + 10;
end
