function results = DropoutBO(fun,vars,varargin)
% Bayesian Optimization Using Dropout
%
% Syntax:
%   results = DropoutBO(fun,vars);
%   results = DropoutBO(fun,vars,'propertyname','propertyvalue',...)
%   
% Description:
%   Bayesian Optimization algorithm using Dropout b< Li et al
%   (High Dimensional Bayesian Optimization Using Dropout, 2018)
%   Each iteration only d varaibles are optimized.
%
% Input:
%   fun:    function handle, as fun = @(x) yourFunction(x.y,x.y,...)
%   vars:   variables as optimizableVariable, e.g [x,y,...] 
%
% Propertyname/-value pairs:
%   AcqFun - name of the acquisition function (string) which should be used
%   (default: EI)
%   PRandom - Probaility [0,..,1] to fill dropout varaibles with random vaules. For
%             PRandom=0 all dropout values are set corresponding to the
%             best found solution. Default 0.
%   d - Number of variables remain after dropout
%
% Output:
%   results
%
% used subfunction: setargs
%
% Date: 08. July, 2019
% Author: Michael Werner

% Default values
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [], 'sampleSize', 1000, 'AcqFun', 'EI', 'PRandom', 0, 'd' 2}; 
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
    % Pick d random indexes that surive dropout
    d_idx = randperm(numVar, d);
    % Generate uniformly distributed sample distribution
    s = zeros(d,params.sampleSize);
    for j=1:params.sampleSize
        for l=1:d
            s(l,j) = rand() * (vars(d_idx(l)).Range(2) - vars(d_idx(l)).Range(1)) ...
                            +  vars(d_idx(l)).Range(1);
        end
    end
    % Determine next evaluation point using GP and an acquisition function
    d_x_next = AcqFun(x(d_idx,params.numSeed + (i-1)),s,y(params.numSeed + (i-1)));
    
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




