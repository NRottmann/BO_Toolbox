function x_next = EI(x,s,y)
% Acquisition Function for BO - Expected Improvement
%
% TODO: Add here moe descriptions

% Get mean and variance from the GP model
[mu,sigma] = GP(x,s,y);
sigma = diag(sigma);
% Determine next evaluation point using an acquisition function
tau = max(y);
tmp = (mu - tau) ./ sigma;
alpha = (mu - tau) .* normcdf(tmp) + sigma .* normpdf(tmp);
[~,ID] = max(alpha);
x_next = s(:,ID);

end

