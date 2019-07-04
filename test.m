clear all
close all
clc

initToolbox

x = optimizableVariable('x',[-100,100]);
y = optimizableVariable('y',[-100,100]);
vars = [x,y];

fun = @(vars) test_fun(vars.x,vars.y);

results = BO(fun,vars);