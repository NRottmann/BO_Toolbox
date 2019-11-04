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

classdef TaskAimingPD < handle
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
        function obj = TaskAimingPD(InitialState,SimulationParameter)
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
            Kp = optimizableVariable('Kp', [0,10]);
%             Kd = optimizableVariable('Kd',[0,0.5]);
            policyParameterDefs = [Kp];%, Kd];
            for i=1:2
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
            Kp = policyVector(1);
    %        Kd = policyVector(2);
%             path = [reshape(policyVector(3:end), 2, []), x_data.x_puck];
            path = [reshape(policyVector(2:end), 2, []), x_data.x_puck];
                       
            TimeSamples = t0:obj.SimulationParameter.dt:obj.SimulationParameter.T;
            N = length(TimeSamples);           
                
            % Follow the trajectory of the mallet
            objectiveVal = 0;
            gamma1 = 1;
            gamma2 = 1;
            e = zeros(2, length(t0:obj.SimulationParameter.dt:t1));
            pos = zeros(2, length(t0:obj.SimulationParameter.dt:t1));
            
            
            pos(:, 1 ) = x_data.x_mallet;
            cp_id = 1;
            tic;
            for i=1:1:N
                x_data = obj.KineticModel.getPositionData();
                % Perfom PD controlling while not all cp are approached
                if cp_id <= size(path, 2)
                    % compute error and derivative
                    e(:, i) = path(:, cp_id) - x_data.x_mallet;
%                     ed = diff(e(:, 1:i)')'/obj.SimulationParameter.dt;
                    % compute PD controller
                    u = Kp * e(:, i);% + Kd * ed(:, end);
                else
                    u = [0;0];
                end
   
                % Do the simulation
                obj.KineticModel = obj.KineticModel.SimulateAction(x_data.x_mallet + u*obj.SimulationParameter.dt,3);
                
                % Cost function
                % Get the position of the puck
                x_data = obj.KineticModel.getPositionData();
                pos(:, i+1) = x_data.x_mallet;
                % Add acceleration of the mallet to the costs
                xdd = diff(pos(:, 1:i+1)', 2)';
                objectiveVal = objectiveVal + gamma1 * norm(xdd(:,end)/obj.SimulationParameter.dt);
                
                % add the distance to the goal
                objectiveVal = objectiveVal + gamma2 * norm(x_data.x_puck - obj.FieldParameters.x_goal);
                
                % If we made a goal, end iteration and reduce cost
                if obj.KineticModel.getGoalData().goal_count(1)
                    objectiveVal = objectiveVal - 200;
                    break;
                end
                
                % Count up and visualize if required
                if visualize
                    obj.Visualizer.UpdateVisualizationPID(path, pos(:, 1:i));
                    timeSpend = toc;
                    pause(obj.SimulationParameter.dt-timeSpend);
                end
                
                % check if a checkpoint of the path is reached
                gameData = obj.KineticModel.getGameData();
                if cp_id <= size(path, 2) && norm(path(:, cp_id) - x_data.x_mallet) < gameData.r_mallet/2
                    cp_id = cp_id + 1;
                end
            end
        end       
    end
end
