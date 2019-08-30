function [x_next, alpha, idx] = UCB(x,s,y, varargin)
% Acquisition Function for BO - Upper Confidence Bound
% Syntax:
%   results = UCB(x,s,y);
%   results = UCB(x,s,y,'propertyname','propertyvalue',...)
%   
% Description:
%   Calculates the upper confidence bound for the samples s using a GP
%   fitted on points x and values y and returns the point in s with the higest UCB.
%
% Input:
%   Input:
%   x: measurement points as a matrix with D x n
%   s: sample points as a matrix with D x l
%   y: measurements at points x with n x 1
% 
%   D: Dimension of the input space
%
% Propertyname/-value pairs:
%   kappa - scaling factor of sigma. Default 1
%
% Output:
%   x_next - point of s with highest UCB
%   alpha - UCB values for all points in s
%   idx - the index of x_next in s
%
% used subfunction: setargs
%
% Date: 07. July, 2019
% Author: Michael Werner

defaultargs = {'kappa', 1};
params = setargs(defaultargs, varargin);
kappa = params.kappa;
% Get mean and variance from the GP model
[mu,sigma] = GP(x,s,y, varargin{:});
sigma = diag(sigma);
% Determine next evaluation point using an acquisition function
a  = mu + kappa * sigma;
[~,ID] = max(a);
x_next = s(:,ID);
if nargout > 1
   alpha = a; 
end
if nargout > 2
   idx = ID;
end
end

