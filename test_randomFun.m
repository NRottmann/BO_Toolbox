clear all
close all
clc

initToolbox

%% Optimize
% The function settings
randomFun = RandomFun();
x = optimizableVariable('x',[-10,10]);
y = optimizableVariable('y',[-10,10]);
vars = [x,y];
fun = @(vars) randomFun.paramFun(vars.x,vars.y);

% Params for the optimization
sampleSize = 1000;
numVar = 2;
numSeed = 2;
maxIter = 30;
iter = 100;

% Initialize storage capacity
resultsHistoryBO = zeros(maxIter+numSeed,iter);
resultsHistoryHIBO = zeros(maxIter+numSeed,iter);
resultsMaxHistoryBO = zeros(maxIter+numSeed,iter);
resultsMaxHistoryHIBO = zeros(maxIter+numSeed,iter);
resultsMaxBO = zeros(iter,1);
resultsMaxHIBO = zeros(iter,1);

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
    resultsBO = BO(fun,vars,'maxIter', maxIter, 'numSeed', numSeed, 'seedPoints', x0, 'sampleSize', sampleSize);
    resultsHistoryBO(:,j) = resultsBO.valueHistory;
    resultsMaxHistoryBO(:,j) = resultsBO.maxValueHistory;
    resultsMaxBO(j) = resultsBO.bestValue;
    
    results = HIBO_forRand(fun,vars,'maxIter', maxIter, 'numSeed', numSeed, 'seedPoints', x0, 'sampleSize', sampleSize);
    resultsHistoryHIBO(:,j) = results.valueHistory;
    resultsMaxHistoryHIBO(:,j) = results.maxValueHistory;
    resultsMaxHIBO(j) = results.bestValue;
      
    disp(j)
end


%% Plots
figure;
plot(mean(resultsHistoryBO,2))
hold on
plot(mean(resultsHistoryHIBO,2))
legend('BO','HIBO')

figure;
plot(mean(resultsMaxHistoryBO,2))
hold on
plot(mean(resultsMaxHistoryHIBO,2))
legend('BO','HIBO')