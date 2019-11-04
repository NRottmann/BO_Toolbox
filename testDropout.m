clear all
close all
clc

initToolbox

%% Go 
f = Schwefel(10);
results = DropoutBO(@f.call,f.vars, 'minimize', f.minimize);
resultsLearning = DropoutBOLearning(@f.call,f.vars, 'minimize', f.minimize);


%% Plot
plot(1:33, results.maxValueHistory)
hold on
plot(1:33, resultsLearning.maxValueHistory)
legend('normal','learning')

