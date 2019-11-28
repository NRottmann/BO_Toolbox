function [K,varargout] = se_kernel_var_cond(X1,X2,varargin)
% Squared-exponential covariance function (Kernel) with equal characteristic length
% in every dimension, it is a measure of similarity
%
% Formular:
%   k(x,x') = sigma_f^2 * exp(-(1/2) * (x - x')^T Delta (x - x')) + sigma_n^2 kroneckerDelta_qp
%   
%           Delta: diagonal matrix of dimension D with diag(l1^2,l2^2,...)
%
% Syntax:
%   [K] = se_kernel(X1,X2)
%   [K,optionalOutputs] = se_kernel(X1,X2,'propertyname','propertyvalue',...)
%
% Description:
%   Generates a Covariance Matrix between the points given in X1 and X2.
%   This covariance matrix can be used for instance for Gaussian Processes
%
% Input: 
%   X1, X2: Matrices of data points with dimension D x n, D x m
%
% Propertyname/-value pairs:
%   CovParam - array of the covariance function parameters
%       2 x 1 array needed:
%       (1): scale factor, sigma_n (default: 1)
%       (2): characteristic diagonal length matrix, Delta (default: [1; 1; ...])
%   struct - gives out the hyperparameters, additional output
%   argument needed (true or false, default: false)
%
% Output:
%   K - covariance matrix (n x m)
%
% Date: 17. August, 2016
% Author: Nils Rottman

% get dimension
D = length(X1(:,1));
vec = ones(1+D,1);

% Default values
defaultargs = {'CovParam', vec,'struct',0, 'num_feature', 0}; 
params = setargs(defaultargs, varargin);

% look if hyperparameters are wanted
if params.struct == true
    varargout{1} = params.CovParam;
end 

% error checking
if iscell(params.CovParam)
    error('covariance parameters are in a cell array')
end
if (length(params.CovParam) ~= (1+D))
    error('covariance parameters have not the right size')
end

% Rewriting parameters in variables
sigma_n = params.CovParam(1);
DELTA = diag((params.CovParam(2+D-params.num_feature:end).^2).^(-1));

% computing the dimension of the kernel matrix
l1 = length(X1(1,:));
l2 = length(X2(1,:));

% providing the kernel matrix for effcient computing
K = zeros(l1,l2);

% Computing the kernel for each element of the kernel matrix
for i=1:1:l1
    for j=1:1:l2
%         if i <= l1-params.num_feature && j <= l2-params.num_feature
%             continue
%         end
        dist = (X1((size(X1,1)-params.num_feature+1):end, i) -...
               X2((size(X2,1)-params.num_feature+1):end, j));
        K(i,j) = sigma_n^2 * exp(-(1/2) * dist' * DELTA * dist);
    end
end
    
end

