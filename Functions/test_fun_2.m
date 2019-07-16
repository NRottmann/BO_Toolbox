function [z] = test_fun_2(x,y)
% 2-dmensional test function for Bayesian Optimization which has only one
% effective dimension

v = [x; y];
T = [0.7071, 0.7071];
val = T * v; 

a = [0.1 1 1];
z = -a(1)*val^2 + a(2)*val + a(3);


end

