%% Clear everything
clear all
close all
clc

%% Generate a random field
B = 10;
stepSize = 0.1;

x = -B:stepSize:B;
y = -B:stepSize:B;

N = length(x);
M = length(y);

[X,Y] = meshgrid(x,y);
Z = gamrnd(1,1,N,M);

%% Plot the field
figure;
surf(X,Y,Z)

%% Generate the feature
numValues = N*M;
x_f = linspace(-10,10,numValues);
z_f = zeros(numValues,1);
X_idx_F = zeros(numValues,2);
Z_tmp = Z;
for i=1:numValues
    [maxCol,yIdx] = max(Z_tmp);
    [maxValue,xIdx] = max(maxCol);
    X_idx_F(i,:) = [X(yIdx(xIdx),xIdx), Y(yIdx(xIdx),xIdx)];
    z_f(i) = Z_tmp(yIdx(xIdx),xIdx);
    Z_tmp(yIdx(xIdx),xIdx) = -inf;
end

figure;
plot(x_f,z_f)

%% Save
save('randomData','X','Y','Z','x_f','z_f','X_idx_F');



