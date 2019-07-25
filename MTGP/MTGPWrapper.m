function [mu,sigma, params] = MTGPWrapper(x, s, y, params)
% MTGPWrapper 
%
% Syntax:
%   [mu,Sigma] = MTGPWrapper(x,s,y);
%   
% Description:
%  
%
% Input:
%   x: measurement points as a matrix with D x n
%   s: unknown points as a matrix with D x l
%   y: measurements at points x with d x n
% With
%   D: Dimension of the input space
%   d: Dimension of the output space (number of task)
%
% Output:
%   mu - (d x l) matrix of mean predictions at points s
%   sigma - (d x l) matrix of MTGP Variance at points s
%
% Date: 25. July, 2019
% Author: Michael Werner

n = size(x, 2);
d = size(y, 1);
covfunc_x = {'covSEard'}; % TODO as parameter
irank=d;
%% prepare data
x_train = x';
x_predict = s';
y_train = y';
y_train = y_train(:);

idx_kf = repmat(1:d, n,1);
idx_kf = idx_kf(:);
idx_kx = repmat((1:n)',d,1);
idx_kx = idx_kx(:);
nx = ones(n*d, 1);

data = {covfunc_x, x_train, y_train, d, irank, nx, idx_kf, idx_kx};

%% learn MTGP

% Hyper-parameter learning
if isempty(params)
    [logtheta_all, deriv_range] = init_mtgp_default(x_train, covfunc_x, d, irank);
    [logtheta_all, ~] = learn_mtgp(logtheta_all, deriv_range, data);
else
    logtheta_all = params;    
end

%% Predict
[mu, sigma] = predict_mtgp_all_tasks(logtheta_all, data, x_predict);

mu = mu';
sigma = sigma';
if nargout >= 3
    params = logtheta_all;
end

