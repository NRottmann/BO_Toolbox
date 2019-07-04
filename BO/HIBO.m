function results = HIBO(fun,vars,varargin)
% Bayesian Optimization
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
%   AcqFun - name of the acquisition function (string) which should be used
%   (default: )
%
% Output:
%   results
%
% used subfunction: setargs
%
% Date: 02. July, 2019
% Author: Nils Rottmann

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'sampleSize', 1000}; 
params = setargs(defaultargs, varargin);

% Get number of variables to optimize
numVar = length(vars);

% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:numSeed
    x_fun = struct();
    for j=1:numVar
        x_fun.(vars(j).Name) = rand() * (vars(j).Range(2) - vars(j).Range(1)) ...
                            +  vars(j).Range(1);
        x(j,i) = x_fun.(vars(j).Name);
    end
    y(i) = fun(x_fun);
end

% We iterate over maxIter iterations
for i=1:params.maxIter
    % Generate uniformly distributed sample distribution
    s = zeros(numVar,params.sampleSize);
    for j=1:params.sampleSize
        for l=1:numVar
            s(l,j) = rand() * (vars(l).Range(2) - vars(l).Range(1)) ...
                            +  vars(l).Range(1);
        end
    end
    % Get mean and variance from the GP model
    [mu,sigma] = GP(x(:,params.numSeed + (i-1)),s,y(params.numSeed + (i-1)));
    sigma = diag(sigma);
    % Determine next evaluation point using an acquisition function
    tau = max(y(1:params.numSeed + (i-1)));
    tmp = (mu - tau) ./ sigma;
    alpha = (mu - tau) .* normcdf(tmp) + sigma .* normpdf(tmp);
    [~,ID] = max(alpha);
    x_next = s(:,ID);
    % Get the next function value
    x_fun = struct();
    for j=1:numVar
        x_fun.(vars(j).Name) = x_next(j);
        x(j,params.numSeed + i) = x_next(j);
    end
    y(params.numSeed + i) = fun(x_fun);
end

% Give back the results
results.valueHistory = y;
[y_max,id_max] = max(y);
results.bestValue = y_max;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end

end




