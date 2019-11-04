% The MIT License (MIT)
% Copyright (c) 2019 Nils Rottmann
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.

% Task Specification for Air Hockey
% 
% Author: Nils Rottmann
% E-Mail: Rottmann@rob.uni-luebeck.de
%
% Date: 12.04.2019

classdef TaskAiming < handle
   	properties (Access = public)
        % Task parameters
        InitialState;       % Initial configuration of the kinematic model
                            % puck, mallet position and velocity and goals
        FieldParameters;    % Goal positions, width and length of the field
        
        % Simulation parameter
        SimulationParameter;    % dt and T
        
        % Public Members
        KineticModel;
        Visualizer;
    end
    
  	methods (Access = public)
        % CONSTRUCTOR
        function obj = TaskAiming(InitialState,SimulationParameter)
            % Generate instance of class
            obj.KineticModel = AirHockeySim.KineticModel();
            obj.Visualizer = AirHockeySim.Visualizer(obj.KineticModel);
            
            % Get game data
            data = obj.KineticModel.getGameData();
            obj.FieldParameters.x_goal = [0; (data.l_field/2 + data.r_puck)];
            obj.FieldParameters.w_field = data.w_field;         
            obj.FieldParameters.l_field = data.l_field;
            obj.FieldParameters.r_puck = data.r_puck;
            obj.FieldParameters.r_mallet = data.r_mallet;
            
            % Store initial state and simulation parameter
            obj.InitialState = InitialState;
            obj.SimulationParameter = SimulationParameter;
            
            % Set state and simulation parameter in the kinematic model
            obj.KineticModel = obj.KineticModel.SetInitialState(InitialState);
            obj.KineticModel = obj.KineticModel.SetTimeStep(SimulationParameter.dt);
        end
        
        % Create the optimization task
        function [policyParameterDefs, objectiveFunc] = CreateOptimizationTask(obj)
            % The function which have to be minimized
            objectiveFunc = @(action)obj.evaluatePolicy(action,0);
            
            % Define the policy parameter which has to be optimized
            dr = 3*obj.FieldParameters.r_mallet;
            dv = 2;
            alpha = optimizableVariable('alpha', [-pi, pi]);
            vx1 = optimizableVariable('vx1',[-dv,dv]);
            vy1 = optimizableVariable('vy1',[0,2*dv]);
            t1 = optimizableVariable('t1',[0,obj.SimulationParameter.T]);
            policyParameterDefs = [alpha,vx1,vy1,t1];
        end
        
        % Returns the cost depending on the given policy
        function [objectiveVal] = evaluatePolicy(obj, policy, visualize)
            % Set model to initial state for simulation evaluation
            obj.KineticModel = obj.KineticModel.SetInitialState(obj.InitialState);
            
            if visualize
                if obj.Visualizer.window == 0
                    obj.Visualizer = obj.Visualizer.InitializeVisualization;
                    pause(0.5)
                end
            end
            
            % Get the policy parameters
            if istable(policy)
                policyVector = table2array(policy);
            else
                policyVector = policy;
            end 
            
            % Allocate initial conditions
            x0 = obj.InitialState.x_mallet(1);
            y0 = obj.InitialState.x_mallet(2);
            vx0 = obj.InitialState.v_mallet(1);
            vy0 = obj.InitialState.v_mallet(2);
            t0 = 0;
            t0_tmp = t0;
            
            % Allocate policy         
            alpha= policyVector(1);
            x1 = obj.InitialState.x_puck(1) + obj.FieldParameters.r_puck * cos(alpha);
            y1 = obj.InitialState.x_puck(2) + obj.FieldParameters.r_puck * sin(alpha);
            vx1 = policyVector(2);
            vy1 = policyVector(3);
            t1 = policyVector(4);
                 
            TimeSamples = t0:obj.SimulationParameter.dt:t1;
            n = length(TimeSamples);
            t = 1;
            q = zeros(2, length(t0:obj.SimulationParameter.dt:t1));
            qd = zeros(2, length(t0:obj.SimulationParameter.dt:t1));
            qdd = zeros(2, length(t0:obj.SimulationParameter.dt:t1));
            while 1
                % Define the trajectory of the mallet
                Waypoints = [x0,x1;y0,y1];
                TimeOfArrival = [t0_tmp; t1];
                Velocities = [vx0,vx1; vy0,vy1];
                TimeSamples = t0_tmp:obj.SimulationParameter.dt:t1;
                [q_tmp,qd_tmp,qdd_tmp,~] = cubicpolytraj(Waypoints,TimeOfArrival,...
                    TimeSamples,'VelocityBoundaryCondition',Velocities);
                q(:, t:end) = q_tmp;
                qd(:, t:end) = qd_tmp;
                qdd(:, t:end) = qdd_tmp;
                
                % check if mallet collides with a wall
                collision_lr = find(abs(q(1, :)) >= ...
                    (obj.FieldParameters.w_field/2 - obj.FieldParameters.r_mallet), 1);
                collision_ud = find(abs(q(2, :)) >= ...
                    (obj.FieldParameters.l_field/2 - obj.FieldParameters.r_mallet), 1);
                
                if isempty(collision_lr) && isempty(collision_ud)
                    % if no collision is detected continue to trajectory
                    % simulation
                    break 
                end
                
                % TODO: below takes the first collision
                % wouldnt it be better to take the last one?
                
                % If collision detected:
                % split at last collsision and compute
                % tr1 = x0y0t0 -> xc yc, tc ,vx=0, vy=0
                % tr2 = xc yc, tc ,vx=0, vy=0 -> xt yt, tt ,vxt, vyt
                % Reuse genTraj
                
                % in case of a collision with a wall, repeat trajectory
                % calculation from where the collision occured.
                t = min([collision_lr, collision_ud])-1;
                x0 = q(1, t);
                y0 = q(2, t);
                vx0 = qd(1, t);
                vy0 = qd(2, t);
                t0_tmp = (t-1)*obj.SimulationParameter.dt;
                if ~isempty(collision_lr) && collision_lr == t+1
                   vx0 = vx0 * -1;
                end
                if ~isempty(collision_ud) && collision_ud == t+1
                   vy0 = vy0 * -1;
                end                           
            end
            TimeSamples = t0:obj.SimulationParameter.dt:t1;
            TimeSamples = [TimeSamples, ...
                (TimeSamples(end)-obj.SimulationParameter.dt):...
                obj.SimulationParameter.dt:obj.SimulationParameter.T];
            N = length(TimeSamples);
            m = N-n;
            q = [q, repmat(q(:,end),1,m)];
            qd = [qd, repmat(qd(:,end),1,m)];
            qdd = [qdd, repmat(qdd(:,end),1,m)];

            % Follow the trajectory of the mallet
            objectiveVal = 0;
            gamma1 = 1;
            gamma2 = 1;
            tic;
            for i=1:1:N
                objectiveVal = 0;
                % Do the simulation
                obj.KineticModel = obj.KineticModel.SimulateAction(q(:,i),3);
                % Add acceleration of the mallet to the costs
                objectiveVal = objectiveVal + gamma1 * norm(qdd(:,i));
                % Get the position of the puck and add costs
                data = obj.KineticModel.getPositionData();
                objectiveVal = objectiveVal + gamma2 * norm(data.x_puck - obj.FieldParameters.x_goal);
                % If we made a goal, end iteration and reduce cost
                if obj.KineticModel.getGoalData().goal_count(1)
                    objectiveVal = objectiveVal - 1000;
                    break;
                end
                % Count up and visualize if required
                if visualize
                    obj.Visualizer.UpdateVisualization_advanced(x1, y1, vx1, vy1, q);
                    timeSpend = toc;
                    pause(obj.SimulationParameter.dt-timeSpend);
                end
            end
        end       
    end
end
