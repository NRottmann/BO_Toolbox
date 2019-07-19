clear all
close all
clc

initToolbox

%% Plot test function
N = 1000;
xyz = zeros(3,N);
for i=1:N
    xyz(1,i) = rand() * 20 - 10;
    xyz(2,i) = rand() * 20 - 10;
    xyz(3,i) = test_fun_2(xyz(1,i),xyz(2,i));
end
figure;
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'*');

%% Optimize
% x = optimizableVariable('x',[-20,20]);
% y = optimizableVariable('y',[-20,20]);
% vars = [x,y];
% fun = @(vars) test_fun_2(vars.x,vars.y);
x1 = optimizableVariable('x1',[-20,20]);
x2 = optimizableVariable('x2',[-20,20]);
x3 = optimizableVariable('x3',[-20,20]);
x4 = optimizableVariable('x4',[-20,20]);
x5 = optimizableVariable('x5',[-20,20]);
x6 = optimizableVariable('x6',[-20,20]);
x7 = optimizableVariable('x7',[-20,20]);
x8 = optimizableVariable('x8',[-20,20]);
x9 = optimizableVariable('x9',[-20,20]);
x10 = optimizableVariable('x10',[-20,20]);
vars = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10];
fun = @(vars) test_fun_3(vars.x1,vars.x2,vars.x3,vars.x4,vars.x5,vars.x6,vars.x7,vars.x8,vars.x9,vars.x10);

numVar = 10;
numSeed = 2;
maxIter = 10;
iter = 20;
resultsHistoryBO = zeros(maxIter+numSeed,iter);
resultsHistoryHIBO = zeros(maxIter+numSeed,iter);
resultsMaxHistoryBO = zeros(maxIter+numSeed,iter);
resultsMaxHistoryHIBO = zeros(maxIter+numSeed,iter);
resultsMaxBO = zeros(iter,1);
resultsMaxHIBO = zeros(iter,1);
resultsFeaturesHIBO = cell(iter,1);
for j=1:1:iter
    % Generate random starting points
    x0 = zeros(numVar,numSeed);
    for ii=1:numSeed
        for jj=1:numVar
            x0(jj,ii) = rand() * (vars(jj).Range(2) - vars(jj).Range(1)) ...
                            +  vars(jj).Range(1);
        end
    end
    
    % Do the algorithm
    resultsBO = BO(fun,vars,'maxIter', maxIter, 'numSeed', numSeed, 'seedPoints', x0);
    resultsHistoryBO(:,j) = resultsBO.valueHistory;
    resultsMaxHistoryBO(:,j) = resultsBO.maxValueHistory;
    resultsMaxBO(j) = resultsBO.bestValue;
    results = HIBO(fun,vars,'maxIter', maxIter, 'numSeed', numSeed, 'seedPoints', x0);
    resultsHistoryHIBO(:,j) = results.valueHistory;
    resultsMaxHistoryHIBO(:,j) = results.maxValueHistory;
    resultsMaxHIBO(j) = results.bestValue;
    resultsFeaturesHIBO{j} = results.nextFeature;
    disp(j)
end

%% Get value of feature
featureValues = zeros(numSeed + maxIter,1);
for j=1:iter
    for i=numSeed+1:1:numSeed + maxIter
        val = (1/sqrt(2)) * resultsFeaturesHIBO{j}(i-numSeed);
        featureValues(i) = featureValues(i) + test_fun_2(val,val);
    end
end
featureValues = featureValues ./ iter;

%% Plots
figure;
plot(mean(resultsHistoryBO,2))
hold on
plot(mean(resultsHistoryHIBO,2))
hold on
plot(featureValues)
legend('BO','HIBO','Feature')

figure;
plot(mean(resultsMaxHistoryBO,2))
hold on
plot(mean(resultsMaxHistoryHIBO,2))
legend('BO','HIBO')