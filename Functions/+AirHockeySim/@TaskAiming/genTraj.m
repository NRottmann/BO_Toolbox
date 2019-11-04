function [T,t] = genTraj(x0,v0,x1,v1,bp,dt)
% Generates a feasible trajectory, if not the zero trajectory will be given
% back
%
% Inputs:
%   x0,x1:  Initial and end position    [2x1]
%   v0,v1:  Initial and end velocity    [2x1]
%   bp:     Boundary points             [2x4]
%   dt:     Timestep
%
% Outputs:
%   T:      Trajectory
%   t:      Time

% check dimension
% TODO!

nb = 4;     % number boundary points (rectangular)
acc_max = 5;

% Check feasibility, its not perfect, still something TODO!
tx = -inf;
bp = [bp, bp(:,1)];
for i=1:1:nb
    b = [(bp(:,i+1)-bp(:,i)), v1] \ (x1 - bp(:,i));
    if b(2) < 0 && b(2) > tx
        tx = b(2);
    end
end
pm = x1 + v1*tx;            % Contact point with boundary
dx = norm(pm - x1);
acc = 0.5 * (norm(v1)^2 / dx);
if (acc > acc_max)
    disp('Trajectory not feasible!')
    T = zeros(2,1);
    t = 0;
    return
end

% Now we generate the trajectory, therefore we assume that we can fully
% accelerate starting with zero velocity
t_acc = norm(v1)/acc_max;
x01 = x1 - 0.5 * acc_max * t_acc^2 * (v1/norm(v1));     % In between point for the trajectory

t_samples_0 = 0:dt:1;
t_samples_1 = 1+dt:dt:1+t_acc;

[T_0,~,~,~] = cubicpolytraj([x0,x01],[0,1],t_samples_0,'VelocityBoundaryCondition',[v0, zeros(2,1)]);
T_1 = x01 + 0.5 * acc_max * (v1/norm(v1)) .* ((t_samples_1-1).^2);

t = [t_samples_0, t_samples_1];
T = [T_0, T_1];

end

