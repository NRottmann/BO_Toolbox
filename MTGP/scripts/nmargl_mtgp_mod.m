function [nl] = nmargl_mtgp_mod(logtheta_all, covfunc_x, x, y,...
				      m, irank, nx, ind_kf, ind_kx)
% Marginal likelihood for multi-task Gaussian Processes
% 
% [nl] = nmargl_mtgp(logtheta_all, covfunc_x, x, y,...
%                	      m, irank, nx, ind_kf, ind_kx)
%
%
% nl = nmargl_mtgp(logtheta_all, ...) Returns the marginal negative log-likelihood
% logtheta_all: Vector of all parameters: [theta_lf; theta_x; sigma_l]
%                - theta_lf: the parameter vector of the
%                   cholesky decomposition of k_f
%                - theta_x: the parameters of K^x
%                - sigma_l: The log of the noise std deviations for each task
% covfunc_x   : Name of covariance function on input space x
% x           : Unique input points on all tasks 
% y           : Vector of target values
% m           : The number of tasks
% irank       : The rank of K^f 
% nx          : number of times each element of y has been observed 
%                usually nx(i)=1 unless the corresponding y is an average
% ind_kx      : Vector containing the indexes of the data-points in x
%                which each observation y corresponds to
% ind_kf      : Vector containing the indexes of the task to which
%                each observation y corresponds to
%

% Author: Edwin V. Bonilla
% Last update: 23/01/2011
% Modified: Michael Werner
% Date 09.08.2019


% *** General settings here ****
config = get_mtgp_config();
MIN_NOISE = config.MIN_NOISE;
% ******************************

if ischar(covfunc_x), covfunc_x = cellstr(covfunc_x); end % convert to cell if needed

D = size(x,2); %  Dimensionality to be used by call to covfunc_x
n = length(y); 

ltheta_x = eval(feval(covfunc_x{:}));     % number of parameters for input covariance

nlf = irank*(2*m - irank +1)/2;        % number of parameters for Lf
vlf = logtheta_all(1:nlf);             % parameters for Lf

theta_lf = vlf; 
Lf = vec2lowtri_inchol(theta_lf,m,irank);

theta_x = logtheta_all(nlf+1:nlf+ltheta_x);                     % cov_x parameters
sigma2n = exp(2*logtheta_all(nlf+ltheta_x+1:end));              % Noise parameters
Sigma2n = diag(sigma2n);                                        % Noise Matrix
Var_nx = diag(1./nx);

Kx = feval(covfunc_x{:}, theta_x, x);
Kf = Lf*Lf';
K = Kf(ind_kf,ind_kf).*Kx(ind_kx,ind_kx);
K = K + ( Sigma2n(ind_kf,ind_kf) .*Var_nx ); 
Sigma_noise = MIN_NOISE*eye(n);
K = K + Sigma_noise;

try
    L = chol(K)';                        % cholesky factorization of the covariance
    
    alpha = solve_chol(L',y);

    % negative log-likelihood
    nl = 0.5*y'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi);
catch
   nl = Inf; 
end
 





 

