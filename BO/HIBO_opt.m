function results = HIBO_opt(fun,vars,varargin)
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
defaultargs = {'maxIter', 30, 'numSeed', 3, 'seedPoints', [], 'sampleSize', 1000,'numFeature',1,'AcqFun', 'EI'}; 
params = setargs(defaultargs, varargin);

AcqFun = str2func(params.AcqFun);

% Check params
if params.numSeed < 2
    params.numSeed = 2;
    disp('HIBO: numSeed has to be at least 2! Set numSeed to 2!')
end

% Get number of variables to optimize
numVar = length(vars);

% Generate storage capacities
x = zeros(numVar,params.maxIter + params.numSeed);
y = zeros(params.maxIter + params.numSeed,1);
y_max = zeros(params.maxIter + params.numSeed,1);

% Start by generating numSeed seedpoints for the BO algorithm
for i=1:params.numSeed
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
            x_fun.(vars(j).Name) = x(j,i);
            x(j,i) = x_fun.(vars(j).Name);
        end
    end
    y(i) = fun(x_fun);
    y_max(i) = max(y);
end

% We iterate over maxIter iterations
T_History = cell(params.maxIter,1);
f_history = zeros(params.numFeature,maxIter);
for i=1:params.maxIter
    % Generate Features
    % T = linearCombination(x,y,params.numFeature);
    
    T = [1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10), ...
                                    1/sqrt(10), 1/sqrt(10), 1/sqrt(10), 1/sqrt(10)];                          
    f = T * x;
    T_History{i} = T;
    
    % Get the range for the learned features by combining them linearly
    rangeFeature = zeros(params.numFeature,2);
    for j=1:params.numFeature
        for l=1:numVar
            rangeFeature(j,1) = rangeFeature(j,1) + vars(l).Range(1) * T(j,l);
            rangeFeature(j,2) = rangeFeature(j,2) + vars(l).Range(2) * T(j,l);
        end
    end
    % Generate uniformly distributed sample distribution for the features
    s_f = zeros(params.numFeature,params.sampleSize);
    for l=1:params.numFeature
        for j=1:params.sampleSize
            s_f(l,j) = rand() * (rangeFeature(l,2) - rangeFeature(l,1)) ...
                            +  rangeFeature(l,1);
        end
    end

    % Determine next feature evaluation point using GP and an acquisition function
    x_next_f = AcqFun(f(:,1:(params.numSeed + (i-1))),s_f,y(1:(params.numSeed + (i-1))));
    f_history(:,i) = x_next_f;
    
    % Generate uniformly distributed sample distribution for the parameters
%     s = zeros(numVar,params.sampleSize);
%     for l=1:numVar
%         for j=1:params.sampleSize
%             s(l,j) = rand() * (vars(l).Range(2) - vars(l).Range(1)) ...
%                             +  vars(l).Range(1);
%         end
%     end
%     % Put both together, to make them dependent
%     x_combined = [x; f];
%     s_combined = [s; x_next_f*ones(1,params.sampleSize)];
%     % Determine next evaluation point using GP and an acquisition function
%     x_combined_next = AcqFun(x_combined(:,1:(params.numSeed + (i-1))),s_combined,y(1:(params.numSeed + (i-1))));
%     x_next = x_combined_next(1:numVar,1);

    x_next = T' * x_next_f;
    
    % TODO: Test, delete later
    % x_next = EI_HIBO(x(:,1:(params.numSeed + (i-1))),s,y(1:(params.numSeed + (i-1))),T);
    
    
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
results.featureHistory = T_History;
results.nextFeature = f_history;
results.valueHistory = y;
results.maxValueHistory = y_max;
[y_max,id_max] = max(y);
results.bestValue = y_max;
for j=1:numVar
    results.bestParams.(vars(j).Name) = x(j,id_max);
end

end

% Feature generation, TODO: Write some documentation and put this in
% own function
function T = linearCombination(x,y,numFeature)
% Feature generation using linear combination

% Get data length and the number of variables
N = length(x(1,:));
numVar = length(x(:,1));

% Get the number of possible combinations and generate matrix
M = 0;
for i=1:(N-1)
    M = M + i;
end
dX = zeros(numVar,M);
dY = zeros(1,M);

% Get the difference
idx = 1;
for i=1:N
    for j=(i+1):N
        dX(:,idx) = x(:,i) - x(:,j);
        dY(:,idx) = y(i) - y(j);
        idx = idx + 1;
    end
end

% use genetic algorithms
options = optimoptions(@ga,'Display','off');
param = ga(@summedDifference,numFeature*numVar,options);

% Generate Transformation Matrix
T = zeros(numFeature,numVar);
for l=1:numFeature
    T(l,:) = param(((l-1)*numVar+1):(numVar*l));
    T(l,:) = T(l,:) / norm(T(l,:));
end

    function L = summedDifference(param)
        T_tmp = zeros(numFeature,numVar);
        for ll=1:numFeature
            T_tmp(ll,:) = param(((ll-1)*numVar+1):(numVar*ll));
            T_tmp(ll,:) = T_tmp(ll,:) / norm(T_tmp(ll,:));
        end
        L = 0;
        for jj=1:N
            L = L + norm(T_tmp*dX(:,jj))/abs(dY(jj));
        end
    end
end




