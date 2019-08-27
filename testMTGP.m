close all
clear all

d=2;
D=3;
covfunc_x       = {'covSEard'};
irank=D;
%% Generate Data
x = (0:0.1:4*pi)'; % (nxd)
x = [x (1:numel(x))'];
s = (0.15:0.1:4*pi)';
s = [s (1:numel(s))'];
n = size(x,1);

y = [sin(x(:,1)); sin(1.5*x(:,1)); cos(x(:,1))];
idx_kf = repmat((1:D),n,1);
idx_kf = idx_kf(:);
idx_kx = repmat((1:n)', D,1);
idx_kx = idx_kx(:);

nx = ones(n*D,1);

%% 2. Assigns cell data for learning and prediction
data  = {covfunc_x, x, y, D, irank, nx, idx_kf, idx_kx};

%% 3. Hyper-parameter learning
[logtheta_all, deriv_range] = init_mtgp_default(x, covfunc_x, D, irank);
[logtheta_all, nl]          = learn_mtgp(logtheta_all, deriv_range, data);

%% 4. Making predictions at all points on all tasks
[ Ypred, Vpred ] = predict_mtgp_all_tasks(logtheta_all, data, s );

%% plot
y = reshape(y, n, D);
figure()
plot(x(:,1), y(:,1), 'r');
hold on
plot(x(:,1), y(:,2), 'b');
plot(x(:,1), y(:,3), 'g');
plot(s(:,1), Ypred(:,1), 'ro')
plot(s(:,1), Ypred(:,2), 'bo')
plot(s(:,1), Ypred(:,3), 'go')
hold off

figure()
plot(x(:,2), y(:,1), 'r');
hold on
plot(x(:,2), y(:,2), 'b');
plot(x(:,2), y(:,3), 'g');
plot(s(:,2), Ypred(:,1), 'ro')
plot(s(:,2), Ypred(:,2), 'bo')
plot(s(:,2), Ypred(:,3), 'go')
hold off

