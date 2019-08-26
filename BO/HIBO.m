function results = HIBO(fun,vars,varargin)
% HiBO: Hierarchical Acquisition Functions for Bayesian Optimization
%
% Syntax:
%   results = HIBO(fun,vars);
%   results = HIBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Extension of the Bayesian Optimization algorithm. With a feature
%   gernerator a feature space is created. By optimizing the acquisition
%   function in feauture space a current optimal feature is selected. The
%   feature is then used to constrain the parameter space while optimizing
%   a second acquisition function.
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
%   FeatureGenerator - defines how the feautre is created 
%                      (default: 'NeuralNet')
%   CovFunc - the covariance/kernel function (default: 'se_kernel_var')
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
%
% used subfunction: setargs, generateSeedPoints, sampleFromRange
%
% Date: 02. July, 2019
% Author: Nils Rottmann
% Date: 15.8.2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [],...
               'sampleSize', 1000, 'AcqFun', 'EI', 'numFeature', 2,...
               'FeatureGenerator', 'NeuralNet',...
               'CovFunc', 'se_kernel_var', 'minimize', false}; 
params = setargs(defaultargs, varargin);

% Check params
if params.numSeed < 2
    params.numSeed = 2;
    disp('HIBO: numSeed has to be at least 2! Set numSeed to 2!')
end

% Get number of variables to optimize
numVar = length(vars);

% extract functions
Cov = str2func(params.CovFunc);
AcqFun = str2func(params.AcqFun);
f_gen = str2func(params.FeatureGenerator);
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
f_history = zeros(params.numFeature,maxIter);
for i=1:params.maxIter
    % Optimize parameters for prediction in feature space
     [covParam, f_genParam] = optimize(x, y, f_gen, Cov);
    
    % Generate Features
    f = f_gen.getfeature(x, f_genParam);
    
    % Generate uniformly distributed sample distribution for the features
    s_f = sampleFromRange(params.numFeature, params.sampleSize,...
                          f_gen.getbounds());
        
    % Determine next feature evaluation point using GP and an acquisition function
    x_next_f = AcqFun(f(:,1:(params.numSeed + (i-1))),  s_f,...
                      y(1:(params.numSeed + (i-1))),...
                      'CovFunc', params.CovFunc, 'CovParam', covParam);
    f_history(:,i) = x_next_f;
    
    % Generate uniformly distributed sample distribution for the parameters
    s = sampleFromRange(numVar, params.sampleSize, vars);
        
    % Put both together, to make them dependent
    x_combined = [x; f];
    s_combined = [s; x_next_f*ones(1,params.sampleSize)];
    
    % Determine next evaluation point using GP and an acquisition function
    x_combined_next = AcqFun(x_combined(:,1:(params.numSeed + (i-1))),...
                             s_combined,y(1:(params.numSeed + (i-1))),...
                             'CovFunc', params.CovFunc);
    x_next = x_combined_next(1:numVar,1);   
    
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
results.nextFeature = f_history;
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

