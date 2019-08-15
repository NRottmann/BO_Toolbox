% This script takes the 2D-Barnin Function and emedded it into a D=25
% dimensional space. The function is then minimized with BO (blue),
% REMBO (red) and DropoutBO (green). Please notize the sign changes due to
% the algorithm maximize!.

clear all
close all
clc

initToolbox

D = 25;
n = 1; % number if 

branin = Branin();
vars = branin.vars;
for i=length(branin.vars):D
   vars = [vars, optimizableVariable(strcat("var", num2str(i)),[0,15])];
end

fun = @(vars) -branin.call(vars.x,vars.y);

N = 200;

valuesREMBO = zeros(N+3,n);
for i=1:n
    results = REMBO(fun,vars, 'numFeatures', 4, 'maxIter', N);
    valuesREMBO(:, i) = results.valueHistory;
end

valuesBO = zeros(N+3,n);
for i=1:n
    results = HIBO(fun,vars,'numFeature',4, 'maxIter', N);
    valuesBO(:, i) = results.valueHistory;
end

valuesHIBO = zeros(N+3,n);
for i=1:n
    results = BO(fun,vars, 'maxIter', N);
    valuesHIBO(:, i) = results.valueHistory;
end

valuesDropout = zeros(N+3,n);
for i=1:n
    results = DropoutBO(fun,vars, 'maxIter', N);
    valuesDropout(:, i) = results.valueHistory;
end

%% Plots
figure;
plot(-mean(valuesBO,2))
hold on
plot(-mean(valuesHIBO,2))
plot(-mean(valuesREMBO,2))
plot(-mean(valuesDropout,2))
legend('BO','HIBO','REMBO','Dropout')
