function results = NNBO(fun,vars,varargin)
% NNBO: Bayesian Optimization with an NN-Kernel
%
% Syntax:
%   results = NNBO(fun,vars);
%   results = NNBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Extension of the Bayesian Optimization algorithm. The parameterspace
%   will be reduced by a neural network. The generated feature space is then 
%   is then used by a GP and a Acquisition function to find a next good
%   feature. the corresponding parameter is then applied to the function
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
%   numFeature - number of features (default: 2)
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
% used subfunction: setargs, generateSeedPoints, sampleFromRange
%
% Date: 30. August, 2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [],...
               'sampleSize', 1000, 'AcqFun', 'EI', 'numFeature', 2,...
               'CovFunc', 'se_kernel_var', 'minimize', false}; 
params = setargs(defaultargs, varargin);

% Check params
if params.numSeed < 2
    params.numSeed = 2;
    disp('MYBO: numSeed has to be at least 2! Set numSeed to 2!')
end

% Get number of variables to optimize
numVar = length(vars);

% extract functions
Cov = str2func(char(params.CovFunc));
AcqFun = str2func(char(params.AcqFun));
%f_gen = str2func('NeuralNet');
f_gen = str2func('LinearFeature');
f_gen = f_gen(numVar, params.numFeature, vars);

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
    x_iter = x(:, 1:(params.numSeed + (i-1)));
    y_iter = y(1:(params.numSeed + (i-1)));
    % Optimize parameters for prediction in feature space
    [covParam, f_genParam] = optimize(x_iter, y_iter, f_gen, Cov);
    
    % Generate Features
    f = f_gen.getfeature(x_iter, f_genParam);
    
    % Generate uniformly distributed sample distribution
    s = sampleFromRange(numVar, params.sampleSize, vars);
    
    % Compute features for the sample points
    s_f = f_gen.getfeature(s, f_genParam);
    
    % Determine next feature evaluation point using GP and an acquisition function
    [~, ~, idx_next] = AcqFun(f,  s_f, y_iter,...
                              'CovFunc', params.CovFunc,...
                              'CovParam', covParam);
    
    % Extract next evaluation point
    x_next = s(:, idx_next);
    
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
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end

end

function [covParam, f_genParam] = optimize(x, y, f_gen, Cov)
% Get number of hyperparameters
f_x = f_gen.getfeature(x, rand(f_gen.num_param,1));
[~,paramTest] = Cov(f_x,f_x,'struct',true);
nvarCov = length(paramTest);

% start optimizaiton
params = optimizeParameter(@logLikelihood, [nvarCov, f_gen.num_param]);
    function L = logLikelihood(param)
        % extract parameters
        covParam = cell2mat(param(1));
        f_gen_param = cell2mat(param(2));
        % genertate features
        f_x = f_gen.getfeature(x, f_gen_param);
        % compute logLikelihood for prediction in feature space
        K = Cov(f_x,f_x,'CovParam',covParam);
        L = (1/2) * (y'/(K + 0.01*eye(length(K(:,1))))) * y...
            + (1/2) * log(norm(K));
    end

% extract parameters
covParam = cell2mat(params(1));
f_genParam = cell2mat(params(2));
end

