function results = DropoutBOLearning(fun,vars,varargin)
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
%   minimize - set true to minimize a function (default: false)
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
               'sampleSize', 1000, 'AcqFun', 'EI', 'PRandom', 0.15,...
               'numFeature' 2, 'CovFunc', 'se_kernel_var',...
               'minimize', false}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(char(params.AcqFun));
PRandom = min(max(params.PRandom, 0), 1);
d = params.numFeature;

% Get number of variables to optimize
numVar = length(vars);

% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
[x, y, y_max] = generateSeedPoints(fun, x, y, y_max, vars,...
                                   params.numSeed, params.seedPoints,...
                                   params.minimize);

% Define importance vector for the dimensions (at the beginning all
% dimensions are equally important)
importance = zeros(numVar,1);
randBase = 1:1:numVar;
                               
% We iterate over maxIter iterations
for i=1:params.maxIter
    % Generate a random vector, for ascending order
    randVec = randBase.^((i-1)/50);
    randVec = randVec/sum(randVec);
    for ii=1:1:numVar-1
        randVec(ii+1) = randVec(ii) + randVec(ii+1);
    end
    % Pick d random indices that survive dropout
    idx_tmp = ones(d,1);
    for ii=1:1:d
        randNum = rand();
        for jj=2:1:numVar
            if randNum > randVec(jj)
                idx_tmp(ii) = jj;
            end
        end
        % Check if idx already used
        used = false;
        for ll=1:1:ii-1
            if(idx_tmp(ii) == idx_tmp(ll))
                used = true;
            end
        end
        if used
            ii = ii - 1;
        end
    end
    % Connect picked indices to importance ordering
    [~,idx_ordering] = sort(importance);
    d_idx  = zeros(d,1);
    for ii=1:1:d
        d_idx(ii) = idx_ordering(idx_tmp(ii));
    end
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
            x_fun.(vars(j).Name) = rand() * (vars(j).Range(2) - ...
                vars(j).Range(1)) + vars(j).Range(1);
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
    if params.minimize
       y(params.numSeed + i) = -y(params.numSeed + i); 
    end
    y_max(params.numSeed + i) = max(y(1:(params.numSeed + i)));
    % Adjust importance vector
    importance(d_idx) = importance(d_idx) + abs(mean(y((params.numSeed + i - 1))) - y(params.numSeed + i));
end

% Give back the results
[ymax,id_max] = max(y);
if params.minimize
    results.valueHistory = -y;
    results.maxValueHistory = -y_max;
    results.bestValue = -ymax;
else
    results.valueHistory = y;
    results.maxValueHistory = y_max;
    results.bestValue = ymax;
end
results.paramHistory = x;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end
end