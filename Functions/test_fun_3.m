function [z] = test_fun_3(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)

v = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];
T = [1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), ...
                                    1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10)];

val = T * v;

a = [0.1 1 1];
z = -a(1)*val^2 + a(2)*val + a(3);

end

