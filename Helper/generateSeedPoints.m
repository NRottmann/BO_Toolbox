function [x, y, y_max] = generateSeedPoints(func, x, y, y_max, vars, numSeed, seedPoints, minimize)
%GENERATESEEDPOINTS generate seed points and evalutes func the those points
%
% Syntax:
%   [x, y, y_max] = generateSeedPoints(func, x, y, y_max, vars, numSeed, seedPoints)
%   
% Description:
%   Generate numSeed seed points based on the given variables vars. 
%    If seed points are given they are copied. The seed points are then
%    applied at func.
%
% Input:
%   func - function used for evaluation
%   x - storage for the seed points
%   y - storage for the function values
%   y_max - maximum over all evaluated function values
%   vars - array of optimizableVariable
%   numSeed - number of seed points
%   seedPoints - given seed points, size: length(var) x numSeed
%   minimize - set true if func has to minimize, else false
%
%
% Output:
%   x - storage for the seed points
%   y - storage for the function values
%   y_max - maximum over all evaluated function values
%
% Date: 12.8.2019
% Author: Michael Werner

numVar = length(vars);
% create numSeed seed points and evaluate objective function at those points
for i=1:numSeed
    x_fun = struct();
    if isempty(seedPoints)
        % if no seed points are given, randomly sample them
        for j=1:numVar
            x_fun.(vars(j).Name) = rand() * (vars(j).Range(2) - vars(j).Range(1)) ...
                            +  vars(j).Range(1);
            x(j,i) = x_fun.(vars(j).Name);
        end
    else
        % if seed points are given, copy them
        if length(seedPoints(1,:)) ~= numSeed || length(seedPoints(:,1)) ~= numVar
            error('Seed Points have wrong size!')
        end
        for j=1:numVar
            x_fun.(vars(j).Name) = seedPoints(j,i);
            x(j,i) = x_fun.(vars(j).Name);
        end
    end
    y(i) = func(x_fun);
    if minimize
       y(i) = -y(i); 
    end
    y_max(i) = max(y(1:i));
end
end