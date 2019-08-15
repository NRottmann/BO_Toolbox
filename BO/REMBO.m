function results = REMBO(fun,vars,varargin)
% Bayesian Optimization with Random Embedding
%
% Syntax:
%   results = REMBO(fun,vars);
%   results = REMBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Bayesian Optimization with Random Embedding algorithm by Wang et al.
%   (Bayesian optimization in high dimensions via random embeddings, 2013)
%   Bayesian optimization is applied to a feature space. The best found
%   feature is transfered into the input space by appling a randomly but
%   fixed Matrix A.
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
%   (default: 'EI')
%   CovFunc - the covariance/kernel function (default: 'se_kernel_var')
%   numFeature - number of features (default: 2)
%   minimize - set true to minimize a function (default: false)
%
% Output:
%   results
%      results.valueHistory - values received from function evaluation
%      results.maxValueHistory - for each iteration the best function value
%                                so far
%      results.paramHistory - all parameters used for evaluation
%      results.nextFeature - for each iteration the created feature
%      results.bestValue - best seen function value
%      results.bestParams - parameters for the best function value
%      results.A - Matrix A which projects from feature into input space
%
% used subfunction: setargs, generateSeedPoints, sampleFromRange
%
% Date: 15.8.2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [],...
               'sampleSize', 1000, 'AcqFun', 'EI',...
               'CovFunc', 'se_kernel_var', 'numFeatures', 2,...
               'minimize', false}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(params.AcqFun);

% Get number of variables in input space
numVar = length(vars);
numFeature = params.numFeatures;

% extract bounds of input space
range = reshape([vars(:).Range], 2, [] )';
xmin = range(:, 1);
xmax = range(:, 2);
clear range

% generate feature variables
fvars = [];
for i=1:numFeature
    fvars = [fvars, optimizableVariable(strcat("f", num2str(i)),...
                                        [-sqrt(numFeature),...
                                        sqrt(numFeature)])];
end
A = rand(numVar, numFeature);
% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
f = zeros(numFeature,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:numSeed
    f_seed = struct();
    for j=1:numFeature
        f_seed.(fvars(j).Name) = rand() * (fvars(j).Range(2) - fvars(j).Range(1)) ...
                            +  fvars(j).Range(1);
        f(j,i) = f_seed.(fvars(j).Name);
    end
    x_seed = A*cell2mat(struct2cell(f_seed));
    % project x_f into input space bounds if nessessary
    x_seed = min(max(x_seed, xmin), xmax);
    x(:, i) = x_seed;
    for j=1:numVar
        x_fun.(vars(j).Name) = x_seed(j);
    end
    
    y(i) = fun(x_fun);
    if params.minimize
       y(i) = -y(i); 
    end
    y_max(i) = max(y(1:i));
end

% We iterate over maxIter iterations
for i=1:params.maxIter
    % Generate uniformly distributed sample distribution
    s = sampleFromRange(params.numFeature, params.sampleSize, fvars);
        
    % Determine next evaluation point using GP and an acquisition function
    f_next = AcqFun(x(:,params.numSeed + (i-1)),s,...
                    y(params.numSeed + (i-1)), 'CovFunc', params.CovFunc);
    
    % Store feature
    for j=1:numFeature
        f(j,params.numSeed + i) = f_next(j);
    end
    
    % Get the next function value
    x_fun = struct();
    x_next = A*f_next;
    
    % project x_next into input space bounds if nessessary
    x_next = min(max(x_next, xmin), xmax);
    x(:, params.numSeed + i) = x_next;
    for j=1:numVar
        x_fun.(vars(j).Name) = x_next(j);        
    end
    
    y(params.numSeed + i) = fun(x_fun);
    if params.minimize
       y(params.numSeed + i) = -y(params.numSeed + i); 
    end
    y_max(params.numSeed + i) = max(y(1:(params.numSeed + i)));
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
results.nextFeature = f;
results.paramHistory = x;
f_best = x(:, id_max);
for j=1:numFeature
    results.bestFeatures.(fvars(j).Name) = f_best(j);
end
x_best = A*f_best;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x_best(j);
end
results.A = A;
end