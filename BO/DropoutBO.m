function results = DropoutBO(fun,vars,varargin)
% Bayesian Optimization Using Dropout
%
% Syntax:
%   results = DropoutBO(fun,vars);
%   results = DropoutBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Bayesian Optimization algorithm using Dropout by Li et al.
%   (High Dimensional Bayesian Optimization Using Dropout, 2018)
%   Each iteration only d varaibles are optimized.
%
% Input:
%   fun:    function handle, as fun = @(x) yourFunction(x.y,x.y,...)
%   vars:   variables as optimizableVariable, e.g [x,y,...] 
%
% Propertyname/-value pairs:
%   maxIter - number of interations performed by BO (default: 30)
%   numSeed - number of points initialy evaluated (default: 3)
%   seedPoints - given seed Points of size (n x numSeed) with n as the
%                number of variables (default: [])
%   sampleSize - number of samples from the acquisition function 
%                (default: 1000)
%   AcqFun - name of the acquisition function (string) which should be used
%            (default: 'EI')
%   numFeature - number of features (default: 1)
%   PRandom - probability of using Dropout-Random mode (default: 0)
%   CovFunc - the covariance/kernel function (default: 'se_kernel_var')
%
% Output:
%   results
%      results.valueHistory - values received from function evaluation
%      results.maxValueHistory - for each iteration the best function value
%                                so far
%      results.paramHistory - all parameters used for evaluation
%      results.bestValue - best seen function value
%      results.bestParams - parameters for the best function value
%
% used subfunction: setargs, generateSeedPoints
%
% Date: 15.8.2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [],...
               'sampleSize', 1000, 'AcqFun', 'EI', 'PRandom', 0,...
               'numFeature' 2, 'CovFunc', 'se_kernel_var'}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(params.AcqFun);
PRandom = min(max(params.PRandom, 0), 1);
d = params.d;

% Get number of variables to optimize
numVar = length(vars);

% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:numSeed
    x_fun = struct();
    if isempty(params.seedPoints)
        for j=1:numVar
            x_fun.(vars(j).Name) = rand() * (vars(j).Range(2) - vars(j).Range(1)) ...
                                +  vars(j).Range(1);
            x(j,i) = x_fun.(vars(j).Name);
        end
    else
        if length(params.seedPoints(1,:)) ~= params.numSeed ||...
           length(params.seedPoints(:,1)) ~= numVar
            error('Seed Points have wrong size!')
        end
        for j=1:numVar
            x_fun.(vars(j).Name) = params.seedPoints(j,i);
            x(j,i) = x_fun.(vars(j).Name);
        end
    end
    y(i) = fun(x_fun);
    y_max(i) = max(y(1:i));
end

% We iterate over maxIter iterations
for i=1:params.maxIter
    % Pick d random indexes that surive dropout
    d_idx = randperm(numVar, d);
    % Generate uniformly distributed sample distribution
    s = sampleFromRange(params.numFeature, params.sampleSize, vars(d_idx));
    
    % Determine next evaluation point using GP and an acquisition function
    d_x_next = AcqFun(x(d_idx,params.numSeed + (i-1)),s,...
                      y(params.numSeed + (i-1)), 'CovFunc', params.CovFunc);
    
    % construct pararmeter vector for next evaluation
    x_fun = struct();
    % Step one: fill with random values or best parameters so far
    [~,id_max] = max(y);
    r = rand();
    for j=1:numVar
        if r<PRandom
            % Pick random value
            x_fun.(vars(j).Name) = rand() * (vars(d_idx(l)).Range(2) - ...
                vars(d_idx(l)).Range(1)) + vars(d_idx(l)).Range(1);
        else
            % Take value from current best evaluation
            x_fun.(vars(j).Name) = x(j, id_max);
        end
        x(j,params.numSeed + i) = x_fun.(vars(j).Name);
    end
    % Step 2: Override varibles from d_idx
    for j=1:d
       x_fun.(vars(d_idx(j)).Name) = d_x_next(j); 
       x(d_idx(j),params.numSeed + i) = d_x_next(j);
    end
    % Get the next function value
    y(params.numSeed + i) = fun(x_fun);
    y_max(params.numSeed + i) = max(y(1:(params.numSeed + i)));
end

% Give back the results
results.valueHistory = y;
results.maxValueHistory = y_max;
results.paramHistory = x;
[y_max,id_max] = max(y);
results.bestValue = y_max;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end

end