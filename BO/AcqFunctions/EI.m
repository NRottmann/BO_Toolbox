function [x_next, alpha, idx] = EI(x,s,y, varargin)
% Acquisition Function for BO - Expected Improvement
% Syntax:
%   results = EI(x,s,y);
%   
% Description:
%   Calculates the expected improvement for the samples s using a GP
%   fitted on points x and values y and returns the point in s with the higest EI.
%
% Input:
%   Input:
%   x: measurement points as a matrix with D x n
%   s: sample points as a matrix with D x l
%   y: measurements at points x with n x 1
% 
%   D: Dimension of the input space
%
% Output:
%   x_next - point of s with highest EI
%   alpha - PI values for all points in s
%   idx - the index of x_next in s
%
% Date: 07. July, 2019
% Author: Michael Werner

% Get mean and variance from the GP model
[mu,sigma] = GP(x,s,y,varargin{:});
sigma = diag(sigma);
% Determine next evaluation point using an acquisition function
tau = max(y);
tmp = (mu - tau) ./ sigma;
a = (mu - tau) .* normcdf(tmp) + sigma .* normpdf(tmp);
[~,ID] = max(a);
x_next = s(:,ID);
if nargout > 1
   alpha = a; 
end
if nargout > 2
   idx = ID;
end
end

