clear all
close all
clc

initToolbox

%% Go 
iter = 100;
numFeatures = 5;
iterBO = 50;
f = Ackley(30);
for i=1:1:iter
    results{i} = DropoutBO(@f.call,f.vars, 'minimize', f.minimize,'numFeature', numFeatures, 'maxIter', iterBO);
    resultsLearning{i} = DropoutBOLearning(@f.call,f.vars, 'minimize', f.minimize,'numFeature', numFeatures, 'maxIter', iterBO);
    disp(i)
end
%% Plot
maxValue = zeros(iterBO+3,iter);
maxValueLearning = zeros(iterBO+3,1);
for i=1:1:iter
    maxValue(:,i) = results{i}.maxValueHistory;
    maxValueLearning(:,i) = resultsLearning{i}.maxValueHistory;
end
mu = mean(maxValue,2);
muLearning = mean(maxValueLearning,2);
sigma = std(maxValue,0,2);
sigmaLearning = std(maxValueLearning,0,2);

f = [(mu+sigma); flip(mu-sigma)];
fLearning = [muLearning+sigmaLearning; flip(muLearning-sigmaLearning)];

figure
fill([1:length(mu), length(mu):-1:1], f, 'blue', 'FaceAlpha', 0.1);
hold on
fill([1:length(mu), length(mu):-1:1], fLearning, 'red', 'FaceAlpha', 0.1);
plot(1:iterBO+3, mu,'b')
plot(1:iterBO+3, muLearning,'r')
legend('normal','learning')

%% Saving the data
% save('dropout_30_5','results','resultsLearning')
