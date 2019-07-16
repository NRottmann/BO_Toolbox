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
defaultargs = {'maxIter', 30, 'numSeed', 3, 'sampleSize', 1000, 'AcqFun', 'EI'}; 
params = setargs(defaultargs, varargin);
AcqFun = str2func(params.AcqFun);
% Get number of variables to optimize
numVar = length(vars);

% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:params.numSeed
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
    for l=1:numVar
        for j=1:params.sampleSize
            s(l,j) = rand() * (vars(l).Range(2) - vars(l).Range(1)) ...
                            +  vars(l).Range(1);
        end
    end
    % Determine next evaluation point using GP and an acquisition function
    x_next = AcqFun(x(:,1:(params.numSeed + (i-1))),s,y(1:(params.numSeed + (i-1))));
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




