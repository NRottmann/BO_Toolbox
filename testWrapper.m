close all
clear all
initToolbox
%% Generate Data
x = (0:0.1:4*pi); 
x = [x; repmat(1:numel(x), 2, 1)];
s = (0.15:0.1:4*pi);
s = [s; repmat(1:numel(s), 2, 1)];

y = [sin(x(1, :)); sin(1.5*x(1, :)); cos(x(1, :)); repmat(1:numel(x(1,:)), 7, 1)];

[mu, sigma] = MTGPWrapper(x, s, y, []);

%% plot
figure()
plot(x(1,:), y(1,:), 'r');
hold on
plot(x(1,:), y(2,:), 'b');
plot(x(1,:), y(3,:), 'g');
plot(s(1,:), mu(1,:), 'ro')
plot(s(1,:), mu(2,:), 'bo')
plot(s(1,:), mu(3,:), 'go')
hold off

figure()
plot(x(2,:), y(1,:), 'r');
hold on
plot(x(2,:), y(2,:), 'b');
plot(x(2,:), y(3,:), 'g');
plot(s(2,:), mu(1,:), 'ro')
plot(s(2,:), mu(2,:), 'bo')
plot(s(2,:), mu(3,:), 'go')
hold off

