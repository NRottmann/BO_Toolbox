function results = ManifoldBO(fun,vars,varargin)
% Manfold Bayesian Optimization
%
% Syntax:
%   results = ManifoldBO(fun,vars);
%   results = ManifoldBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Manifold Bayesian Optimization algorithm
%   The algorithm jointly learns a feature mapping, a GP to model the
%   function and a reconstruction of the input using a manifold 
%   MultiOuputGP. 
%   For futher details pleas refer to:
%   "High-dimensional Bayesian optimization using 
%   low-dimensional feature spaces" from Moriconi et al.
%
% Input:
%   fun:    function handle, as fun = @(x) yourFunction(x.y,x.y,...)
%   vars:   variables as optimizableVariable, e.g [x,y,...] 
%
% Propertyname/-value pairs:
%   AcqFun - name of the acquisition function (string) which should be used
%   (default: )
    % TODO parameters
%
% Output:
%   results
%
% used subfunction: setargs
%
% Date: 24. July, 2019
% Author: Franz Johannes Michael Werner

defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [], 'sampleSize', 1000,...
    'num_hidden', 20, 'feature_size', 2, 'AcqFun', 'EI', 'CovFunc', 'se_kernel_var'}; 
params = setargs(defaultargs, varargin);

% Defining the call for the covariance function
Cov = str2func(params.CovFunc);
AcqFun = str2func(params.AcqFun);

% Create neural net
architecture = [length(vars), params.num_hidden, params.feature_size]; 

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
        if length(params.seedPoints(1,:)) ~= params.numSeed || length(params.seedPoints(:,1)) ~= numVar
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
    x_i = x(:,1:params.numSeed + (i-1));
    y_i = y(1:params.numSeed + (i-1));
    
    % Optimize hyperparameters
    [weights, covParam, mtgpParam] = optimize(x_i, y_i, architecture, Cov);
    
    % Generate uniformly distributed sample distribution
    s = zeros(numVar,params.sampleSize);
    for j=1:params.sampleSize
        for l=1:numVar
            s(l,j) = rand() * (vars(l).Range(2) - vars(l).Range(1)) ...
                            +  vars(l).Range(1);
        end
    end
    f_x = NeuralNet(x_i, weights, architecture);
    f_s = NeuralNet(s, weights, architecture);

    % Determine next evaluation point using GP and an acquisition function
    f_x_next = AcqFun(f_x, f_s, y_i, 'CovParam', covParam);
    
    % map feature into input space
    [x_next, ~] = MTGPWrapper(f_x, f_x_next, x_i, mtgpParam);  
    
    % Get the next function value
    x_fun = struct();
    for j=1:numVar
        x_fun.(vars(j).Name) = x_next(j);
        x(j,params.numSeed + i) = x_next(j);
    end
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

function [netParam, covParam, mtgpParam] = optimize(x,y, architecture, Cov)
    % Optimize parameter for a defined CovFunc, neual net and MTGP

    % Get number of hyperparameters
    nvarNet = NumberOfWeights(architecture);
    f_x = NeuralNet(x, rand(1,nvarNet), architecture);
    [~,paramTest] = Cov(f_x,f_x,'struct',true);
    nvarCov = length(paramTest);
    
    nvar = nvarCov + nvarNet;

    % Use fminunc
    param_0 = 0.1*ones(nvar,1);
    options = optimoptions(@fminunc,'Display','off');
    param = fminunc(@logLikelihood,param_0,options);
        function L = logLikelihood(param)
            covParam = param(1: nvarCov);
            weights = param(nvarCov + 1: end);
            f_x = NeuralNet(x, weights', architecture);
            K = Cov(f_x,f_x,'CovParam',covParam);
            L = (1/2) * y' * inv(K + 0.001*eye(length(K(:,1)))) * y + (1/2) * log(norm(K));
        end
    
    % extract parameters
    covParam = param(1:nvarCov);
    netParam = (param(nvarCov+1 : end))';  
    
    % get parameters for MTGP
    f_x = NeuralNet(x, netParam, architecture);
    [~, ~, mtgpParam] = MTGPWrapper(f_x, f_x, x, []);
end

function [f] = NeuralNet(x, weights, architecture, varargin)
% Neual Net
%
% Syntax:
%   results = ManifoldBO(x, weights, architecture, varargin);
%   
% Description:
%   Implements the forward propagation of a Neural Net
%
% Input:
%   x: input data with size n x architecture(1)
%   weights: weigths of the neural net 
%   architecture: defines the architecture of the net.
%                 architecture(1) is the input size.
%                 architecture(end) is the output size. 
%
% Output:
%   f
%
% Date: 04. July, 2019
% Author: Franz Johannes Michael Werner

if(length(architecture) < 2)
   error("Need at least 2 layers")
end
if(size(x, 1) ~= architecture(1))
    error("Input size does not match input layer size")
end

y = x';
w = weights;
for layer_idx=1:length(architecture)-1
    W = w(1:architecture(layer_idx) * architecture(layer_idx+1));
    W = reshape(W, architecture(layer_idx), architecture(layer_idx+1));
    b = w(architecture(layer_idx) * architecture(layer_idx+1) + 1 : ...
        (architecture(layer_idx) + 1) * architecture(layer_idx+1));
    a =  y * W + b;
    y =  1.0 ./ (1 + exp(-a));
    
    w = w((architecture(layer_idx) + 1) * architecture(layer_idx+1):end); 
end
f = y';
end

function n_weights = NumberOfWeights(architecture)
% Computes the number of weights for a given architecture.
    n_weights = 0;
    for layer_idx=1:length(architecture)-1
        n_weights = n_weights + (architecture(layer_idx) + 1) * architecture(layer_idx +1);
    end
end