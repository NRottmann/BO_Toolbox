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
%   AcqFun - name of the acquisition function (string) which should be used
%   (default: EI)
%   numFeatures - number of dimensions in feature space
%
% Output:
%   results
%
% used subfunction: setargs
%
% Date: 08. July, 2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [], 'sampleSize', 1000,...
    'numFeatures', 1, 'AcqFun', 'EI'}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(params.AcqFun);

% Get number of variables in input space
numVar = length(vars);
numFeature = params.numFeatures;
% generate feature variables
fvars = [];
for i=1:numFeature
    fvars = [fvars, optimizableVariable(strcat("f", num2str(i)),...
        [-sqrt(numFeature), sqrt(numFeature)])];
end
A = rand(numVar, numFeature);
% Generate storage capacities
x = zeros(numFeature,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:numSeed
    f = struct();
    for j=1:numFeature
        f.(fvars(j).Name) = rand() * (fvars(j).Range(2) - fvars(j).Range(1)) ...
                            +  fvars(j).Range(1);
        x(j,i) = f.(fvars(j).Name);
    end
    x_f = A*cell2mat(struct2cell(f));
    for j=1:numVar
        x_fun.(vars(j).Name) = x_f(j);
    end
    
    y(i) = fun(x_fun);
    y_max(i) = max(y(1:i));
end

% We iterate over maxIter iterations
for i=1:params.maxIter
    % Generate uniformly distributed sample distribution
    s = zeros(numFeature,params.sampleSize);
    for j=1:params.sampleSize
        for l=1:numFeature
            s(l,j) = rand() * (fvars(l).Range(2) - fvars(l).Range(1)) ...
                            +  fvars(l).Range(1);
        end
    end
    % Determine next evaluation point using GP and an acquisition function
    f_next = AcqFun(x(:,params.numSeed + (i-1)),s,y(params.numSeed + (i-1)));
    % Store feature
    for j=1:numFeature
        x(j,params.numSeed + i) = f_next(j);
    end
    % Get the next function value
    x_fun = struct();
    x_next = A*f_next;
    for j=1:numVar
        x_fun.(vars(j).Name) = x_next(j);        
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
f_best = x(:, id_max);
for j=1:numFeature
    results.bestFeatures.(fvars(j).Name) = f_best(j);
end
x_best = A*f_best;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x_best(j);
end

end




