% The MIT License (MIT)
% Copyright (c) 2019 Michael Werner
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
% Author: Michael Werner
% E-Mail: franz.werner@student.uni-luebeck.de
%
% Date: 14.10.2019

classdef TaskAimingSpline < handle
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
        function obj = TaskAimingSpline(InitialState,SimulationParameter)
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
            objectiveFunc = @(policy)obj.evaluatePolicy(policy, 0);
            
            % Define the policy parameter which has to be optimized
            x_data = obj.KineticModel.getGameData();
            gameData = obj.KineticModel.getGameData();
            %dt = optimizableVariable('dt', [0.01,0.05]);
            %policyParameterDefs = dt;
            policyParameterDefs = [];
            for i=1:10
               policyParameterDefs = [policyParameterDefs,...
                   optimizableVariable(strcat('x', num2str(i)), [-x_data.w_field/2 + gameData.r_mallet, x_data.w_field/2 - gameData.r_mallet]),...
                   optimizableVariable(strcat('y', num2str(i)), [-x_data.l_field/2 + gameData.r_mallet, 0])];
            end            
        end
        
        % Returns the cost depending on the given policy
        function [objectiveVal] = evaluatePolicy(obj, policy, visualize)
            % Set model to initial state for simulation evaluation
            obj.KineticModel = obj.KineticModel.SetInitialState(obj.InitialState);
            
            if visualize
                if isempty(obj.Visualizer.window)
                    obj.Visualizer = obj.Visualizer.InitializeVisualization();
                    pause(0.5)
                end
            end
            
            % Get the policy parameters
            if istable(policy)
                policyVector = table2array(policy);
            elseif iscell(policy)
                policyVector = cell2mat(policy);
            else
                policyVector = policy;
            end 
            
            
            x_data = obj.KineticModel.getPositionData();
            % Allocate initial conditions
            t0 = 0;
            t1 = 5;
            % Allocate policy
            %dt = policyVector(2);
            %waypoints = [x_data.x_mallet, reshape(policyVector(2:end), 2, []), x_data.x_puck];
            waypoints = [x_data.x_mallet, reshape(policyVector, 2, []), x_data.x_puck];
            tt = 0:0.01:2+0.01;
            t = 0:(2/11):2;
            xx = spline(t, waypoints(1,:), tt);
            yy = spline(t, waypoints(2,:), tt);

                       
            TimeSamples = t0:obj.SimulationParameter.dt:obj.SimulationParameter.T;
            N = length(TimeSamples);           
                
            % Follow the trajectory of the mallet
            objectiveVal = 0;
            gamma1 = 1;
            gamma2 = 1;            
            tic;
            for i=1:1:N
                x_data = obj.KineticModel.getPositionData();
                if i <= length(xx)
                    new_x = [xx(i); yy(i)]; 
                else
                    new_x = x_data.x_mallet;
                end
   
                % Do the simulation
                obj.KineticModel = obj.KineticModel.SimulateAction(new_x,3);
                
                % Cost function
                % Get the position of the puck
                x_data = obj.KineticModel.getPositionData();
                % Add acceleration of the mallet to the costs
                if i < length(xx)
                    xdd = diff([xx(:, 1:i+1); yy(:, 1:i+1)]', 2)';
                    objectiveVal = objectiveVal + gamma1 * norm(xdd(:,end)/obj.SimulationParameter.dt);
                end
                
                % add the distance to the goal
                objectiveVal = objectiveVal + gamma2 * norm(x_data.x_puck - obj.FieldParameters.x_goal);
                
                % If we made a goal, end iteration and reduce cost
                if obj.KineticModel.getGoalData().goal_count(1)
                    objectiveVal = objectiveVal - 500;
                    break;
                end
                
                % Count up and visualize if required
                if visualize
                    obj.Visualizer.UpdateVisualizationPID(waypoints, [xx; yy]);
                    timeSpend = toc;
                    pause(obj.SimulationParameter.dt-timeSpend);
                end
            end
        end       
    end
end
