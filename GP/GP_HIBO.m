function [mu,sigma] = GP_HIBO(x,s,y,varargin)
% Gaussian Process implemented with zero mean
%
% Syntax:
%   [mu,Sigma] = GP(x,s,y);
%   [mu,Sigma] = GP(x,s,y,'propertyname','propertyvalue',...)
%   
% Description:
%   Gaussian Process implemented with zero mean. Gives a stochastical
%   estimation for the points s for the searched value.
%
% Input:
%   x: measurement points as a matrix with D x n
%   s: unknown points as a matrix with D x l
%   y: measurements at points x with n x 1
% 
%   D: Dimension of the input space
%
% Propertyname/-value pairs:
%   noise - noise of the measurements (default: noise = 0.1)
%   posDef - Parameter for ensuring positive definitness (default: gamma =
%   0.1)
%   CovFunc - name of the covariance function (string) which should be used
%   (default: se_kernel)
%   CovParam - Array of the Covariance Parameters, for further
%   information see Covariance Function description
%
% Output:
%   mu - mean values for points s
%   Sigma - Variance and Covariance regarding the points s
%   time - needed comptutational time
%
% used subfunction: setargs
%
% Date: 01. July, 2019
% Author: Nils Rottmann

% Default values
defaultargs = {'noise', 0.1, 'posDef', 0.1, 'CovFunc', 'se_kernel_var_HIBO', 'CovParam', []}; 
params = setargs(defaultargs, varargin);

% error checking
if size(x,1) ~= size(s,1)
    error('x and s must have the same number of rows')
end

if size(x,2) ~= size(y,1)
    error('The number of columns of x has to be the same as the numbe of rows in y')
end

% check for nans and delete them 
II = any(isnan(x),1);
III = isnan(y);
x(:,II) = [];
y(II)   = [];
x(:,III) = [];
y(III)   = [];

% Defining the call for the covariance function
Cov = str2func(params.CovFunc);

% rewriting parameter
sigma_n = params.noise;

% Creating covariance matrices
if isempty(params.CovParam)  % If no parameters are given take default values
    covParam = covParamOptim(x,y,params.CovFunc);
    [K] = Cov(x,x,'CovParam',covParam);
    [K_s] = Cov(x,s,'CovParam',covParam);
    [K_ss] = Cov(s,s,'CovParam',covParam);
else
    [K] = Cov(x,x,'CovParam',params.CovParam);
    [K_s] = Cov(x,s,'CovParam',params.CovParam);
    [K_ss] = Cov(s,s,'CovParam',params.CovParam);
end

% Ensure positive definiteness
K = 0.1*eye(length(K)) + K;                    

% Calculating mean vector and covariance matrix
I = eye(length(K));
L = chol((K + I*sigma_n^2),'lower');    % Cholesky Tranformation, L satisfies L*L' = K
alpha = (L')\(L\y);                     % avoiding inverse transformation for computational effectiveness
mu = (K_s')*alpha;
beta = L \ K_s;
sigma = K_ss - (beta')*beta;

end

% Log Likelihood Function, TODO: Write some documentation and put this in
% own function
function param = covParamOptim(x,y,CovFunc)
% Optimize parameter for a defined CovFunc

% define handle for covariance function
Cov = str2func(CovFunc);

% Get number of hyperparameters
[~,paramTest] = Cov(x,x,'struct',true);
nvar = length(paramTest);

% Use fminunc
param_0 = ones(nvar,1);
options = optimoptions(@fminunc,'Display','off');
param = fminunc(@logLikelihood,param_0,options);

% % use genetic algorithms
% options = optimoptions(@ga,'Display','off');
% param = ga(@logLikelihood,nvar,options);

    function L = logLikelihood(param)
        K = Cov(x,x,'CovParam',param);
        L = (1/2) * y' * inv(K + 0.001*eye(length(K(:,1)))) * y + (1/2) * log(norm(K));
    end
end

