clear all
close all
clc

initToolbox

x = optimizableVariable('x',[-10,10]);
y = optimizableVariable('y',[-10,10]);
vars = [x,y];

fun = @(vars) test_fun(vars.x,vars.y);

resultsBO = BO(fun,vars);
resultsHIBO = HIBO(fun,vars);


%% Plots
figure;
plot(resultsBO.valueHistory)
hold on
plot(resultsHIBO.valueHistory)
legend('BO','HIBO')

% figure;
% plot(resultsBO.paramHistory(1,:),resultsBO.paramHistory(2,:),'*')