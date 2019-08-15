function results = BO(fun,vars,varargin)
% Bayesian Optimization for maximum search
%
% Syntax:
%   results = BO(fun,vars);
%   results = BO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Bayesian Optimization algorithm
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
%   minimize - set true to minimize a function (default: false)
%   numFeature - only added for compatibility. Not used !!
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
% used subfunction: setargs, generateSeedPoints, sampleFromRange
%
% Date: 02. July, 2019
% Author: Nils Rottmann
% Date: 15.8.2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [],...
               'sampleSize', 1000, 'AcqFun', 'EI',...
               'CovFunc', 'se_kernel_var', 'numFeature', 0,...
               'minimize', false}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(params.AcqFun);

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

% We iterate over maxIter iterations
for i=1:params.maxIter
    % Generate uniformly distributed sample distribution
    s = sampleFromRange(numVar, params.sampleSize, vars);
       
    % Determine next evaluation point using GP and an acquisition function
    x_next = AcqFun(x(:,1:(params.numSeed + (i-1))),s,...
                    y(1:(params.numSeed + (i-1))),...
                    'CovFunc', params.CovFunc);
    
    % Get the next function value
    x_fun = struct();
    for j=1:numVar
        x_fun.(vars(j).Name) = x_next(j);
        x(j,params.numSeed + i) = x_next(j);
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
results.paramHistory = x;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end

end




