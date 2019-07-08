% This script takes the 2D-Barnin Function and emedded it into a D=25
% dimensional space. The function is then minimized with BO (blue),
% REMBO (red) and DropoutBO (green). Please notize the sign changes due to
% the algorithm maximize!.

clear all
close all
clc

initToolbox

D = 25;
n = 10; % number if 

x = optimizableVariable('x',[-5,10]);
y = optimizableVariable('y',[0,15]);
vars = [x,y];
for i=3:D
   vars = [vars, optimizableVariable(strcat("var", num2str(i)),[0,15])];
end

fun = @(vars) -Branin(vars.x,vars.y);
values = zeros(33,n);
for i=1:n
    results = REMBO(fun,vars, 'numFeatures', 4);
    values(:, i) = results.valueHistory;
end
m = mean(values, 2);
s = 0.25*1.96*sqrt(var(values,0,2));
figure(1)
set(gca, 'YScale', 'log')
hold on
plot(1:length(m), -m, 'r', 'DisplayName', 'REMBO');
fill([1:length(m), fliplr(1:length(m))], -[(m+s)', fliplr((m-s)')], 'r', 'FaceAlpha', 0.1);
hold off

values = zeros(33,n);
for i=1:n
    results = BO(fun,vars);
    values(:, i) = results.valueHistory;
end
m = mean(values, 2);
s = var(values,0,2);
hold on
plot(1:length(m), -m, 'b', 'DisplayName', 'BO');
fill([1:length(m), fliplr(1:length(m))], -[(m+s)', fliplr((m-s)')], 'b', 'FaceAlpha', 0.1);
hold off

values = zeros(33,n);
for i=1:n
    results = DropoutBO(fun,vars);
    values(:, i) = results.valueHistory;
end
m = mean(values, 2);
s = var(values,0,2);
hold on
plot(1:length(m), -m, 'g', 'DisplayName', 'BO');
fill([1:length(m), fliplr(1:length(m))], -[(m+s)', fliplr((m-s)')], 'g', 'FaceAlpha', 0.1);
hold off