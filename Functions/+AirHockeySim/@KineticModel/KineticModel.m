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

% Kinetic Model for Air Hockey
% 
% Author: Nils Rottmann
% E-Mail: Rottmann@rob.uni-luebeck.de
%
% Date: 18.03.2019

classdef KineticModel < handle
   	properties (Access = public)
    end
    properties (Access = private)
        x_puck;                     % Positions
        x_mallet;
        
        v_puck;                     % Velocities
        v_mallet;
        
        dt;                         % Simulation time step
        
        r_puck = 0.05;              % Radius puck and mallet
        r_mallet = 0.08;
        
        m_puck = 0.1;               % Masses of the puck and the mallet
        m_mallet = 0.3;
        
        w_field = 1;                % Field width and length
        l_field = 2;
        x_post;                     % Positions of the posts of the goals
        
        w_goal = 0.20;              % WIdth of the goals
        
        g = 9.81;                   % gravitational acceleration
        mu_fr = 0.1;                % friction coefficient
        
        goal_count;                 % Goal counts
    end
    
  	methods (Access = public)
        % CONSTRUCTOR
        function obj = KineticModel()
            % Get initial values
            out = config('AirHockeyInitialSetting');
            obj.x_puck =  out.x_puck;
            obj.x_mallet =  out.x_mallet;
            obj.v_puck =  out.v_puck;
            obj.v_mallet =  out.v_mallet;
            obj.dt =  out.dt;
            obj.goal_count =  out.goal_count;
            
            % Define post positions
            obj.x_post = zeros(2,4);
            count = 1;
            for x_idx = -1:2:1
                for y_idx = -1:2:1
                    obj.x_post(:,count) = [x_idx * obj.w_goal/2; y_idx * obj.l_field/2];
                    count = count + 1;
                end
            end
        end
        
        % Function to simulate the air hockey table
        function [obj] = SimulateAction(obj,u,mode)
            % mode:
            %   1:  Add u (2x1) as a force to the mallet
            %   2:  Add u (2x1) as a force to the puck
            %   3:  Set position of the mallet u (2x1)
            
            % Calculate accelerations
            R = -obj.mu_fr * obj.g;
            phi_puck = atan2(obj.v_puck(2),obj.v_puck(1));
            phi_mallet = atan2(obj.v_mallet(2),obj.v_mallet(1));
            if mode == 1 
                % Generate Accelerations
                a_puck(1,1) = cos(phi_puck) * R;
                a_puck(2,1) = sin(phi_puck) * R;
                a_mallet(1,1) = cos(phi_mallet) * R + (u(1)/obj.m_mallet); 
                a_mallet(2,1) = sin(phi_mallet) * R + (u(2)/obj.m_mallet);
                % Update velocities
                obj.v_puck = obj.v_puck + a_puck * obj.dt;
                obj.v_mallet = obj.v_mallet + a_mallet * obj.dt;
                % Update positions
                obj.x_puck = obj.x_puck + obj.v_puck * obj.dt;
                obj.x_mallet = obj.x_mallet + obj.v_mallet * obj.dt;
            elseif mode == 2
                % Generate Accelerations
                a_puck(1,1) = cos(phi_puck) * R + (u(1)/obj.m_puck);
                a_puck(2,1) = sin(phi_puck) * R + (u(2)/obj.m_puck);
                a_mallet(1,1) = cos(phi_mallet) * R; 
                a_mallet(2,1) = sin(phi_mallet) * R;
                % Update velocities
                obj.v_puck = obj.v_puck + a_puck * obj.dt;
                obj.v_mallet = obj.v_mallet + a_mallet * obj.dt;
                % Update positions
                obj.x_puck = obj.x_puck + obj.v_puck * obj.dt;
                obj.x_mallet = obj.x_mallet + obj.v_mallet * obj.dt;
            elseif mode == 3
                % Calculate acceleration of the puck
                a_puck(1,1) = cos(phi_puck) * R;
                a_puck(2,1) = sin(phi_puck) * R;
                % Update velocity of the puck
                obj.v_puck = obj.v_puck + a_puck * obj.dt;
                % Update position of the puck
                obj.x_puck = obj.x_puck + obj.v_puck * obj.dt;
                % Calculate velocity of the mallet
                obj.v_mallet = (u - obj.x_mallet)/obj.dt;
                % Set position of the mallet
                obj.x_mallet = u;
            else
                error('KineticModel.SimulateAction: Wrong mode chosen!')
            end
            
            % Checks for collisions and restrictions
            % Restrict positions
            if obj.x_mallet(2) > -obj.r_mallet
                obj.x_mallet(2) = -obj.r_mallet;
                obj.v_mallet(2) = 0;
            end
            % Check hitting the side walls
            if abs(obj.x_puck(1)) >= (obj.w_field/2 - obj.r_puck)
                obj.x_puck(1) = sign(obj.x_puck(1)) * (obj.w_field/2 - obj.r_puck);
                obj.v_puck(1) = -obj.v_puck(1);
            end
            if abs(obj.x_mallet(1)) >= (obj.w_field/2 - obj.r_mallet)
                obj.x_mallet(1) = sign(obj.x_mallet(1)) * (obj.w_field/2 - obj.r_mallet);
                obj.v_mallet(1) = -obj.v_mallet(1);
            end
            % Check when the puck approaches the top and bottom wall
            if abs(obj.x_puck(2)) >= (obj.l_field/2 - obj.r_puck)
                % Puck is not between the two goal posts
                if abs(obj.x_puck(1)) > (obj.w_goal/2)
                    obj.x_puck(2) = sign(obj.x_puck(2)) * (obj.l_field/2 - obj.r_puck);
                    obj.v_puck(2) = -obj.v_puck(2);
                % Puck is between the two goal posts
                else
                    count = 1;
                    d = zeros(4,1);
                    for x_idx = -1:2:1
                        for y_idx = -1:2:1
                            d(count) = sqrt((obj.x_post(1,count) - obj.x_puck(1))^2 ...
                                                + (obj.x_post(2,count) - obj.x_puck(2))^2);
                            count = count + 1;
                        end
                    end
                    [dmin,idx] = min(d);
                    if dmin <= obj.r_puck
                        [~,~,obj.v_puck,~] = obj.bounce(obj.v_puck, ...
                            [0; 0],obj.x_puck,obj.x_post(:,idx), ...
                            obj.m_puck,10000,obj.r_puck,obj.r_mallet);
                    end
                end
            end
            % Check when the mallet approaches the top and bottom wall
            if abs(obj.x_mallet(2)) >= (obj.l_field/2 - obj.r_mallet)
                obj.x_mallet(2) = sign(obj.x_mallet(2)) * (obj.l_field/2 - obj.r_mallet);
                obj.v_mallet(2) = -obj.v_mallet(2);
            end
            % Check if the puck hits the mallet
            dx = obj.x_puck - obj.x_mallet;
            dx_norm = norm(dx);    
            if dx_norm < (obj.r_puck + obj.r_mallet)
                % Calculate bounce and set positions
                [obj.x_puck,obj.x_mallet,obj.v_puck,obj.v_mallet] = ...
                    obj.bounce(obj.v_puck,obj.v_mallet,obj.x_puck, ...
                    obj.x_mallet,obj.m_puck,obj.m_mallet, ...
                    obj.r_puck,obj.r_mallet);           
            end
            %check to see if one has scored a goal
            if obj.x_puck(2) > (obj.l_field/2 + obj.r_puck)
                obj.goal_count(1) = obj.goal_count(1) + 1;
                obj.x_puck = [0; 0];
                obj.v_puck = [0; 0];
            elseif obj.x_puck(2)< (-obj.l_field/2 - obj.r_puck)
                obj.goal_count(2) = obj.goal_count(2) + 1;
                obj.x_puck = [0; 0];
                obj.v_puck = [0; 0];
            end 
        end
        
        % Function to access the game data
        function data = getGameData(obj)
            % Function to access the game data for visualization, such as
            % field width or length
            data.r_puck = obj.r_puck;
            data.r_mallet = obj.r_mallet;
            data.w_field = obj.w_field;         
            data.l_field = obj.l_field;
            data.w_goal = obj.w_goal;   
        end
        
        % Function to access the position data of the puck and mallet
        function data = getPositionData(obj)
            % Function to access the position data of the puck and mallet
            data.x_puck = obj.x_puck;
            data.x_mallet = obj.x_mallet; 
        end
        
        % Function to access the goal data
        function data = getGoalData(obj)
            % Function to access the goal data
            data.goal_count = obj.goal_count;
        end
        
        function obj = SetPuckPosition(obj,x,y)
            obj.x_puck = [x; y];
        end
        function obj = SetMalletPosition(obj,x,y)
            obj.x_mallet = [x; y];
        end
        
        function obj = SetToDefault(obj)
            % Get initial values
            out = config('AirHockeyInitialSetting');
            obj.x_puck =  out.x_puck;
            obj.x_mallet =  out.x_mallet;
            obj.v_puck =  out.v_puck;
            obj.v_mallet =  out.v_mallet;
            obj.dt =  out.dt;
            obj.goal_count =  out.goal_count;
        end
        
        function obj = SetInitialState(obj,InitialState)
            % Set initial state
            obj.x_puck =  InitialState.x_puck;
            obj.x_mallet =  InitialState.x_mallet;
            obj.v_puck =  InitialState.v_puck;
            obj.v_mallet =  InitialState.v_mallet;
            obj.goal_count =  InitialState.goal_count;
        end
        
        function obj = SetTimeStep(obj,dt)
            obj.dt = dt;
        end
    end
    
    methods (Static)
        % Function which evaluates the bouncing of the puck
        function [x_I,x_II,v1_I,v1_II] = bounce(v0_I,v0_II,x_I,x_II,m_I,m_II,r_I,r_II)
            % Calculates the velocity of te puck after contact
            % 
            % Assumptions for impact:
            %   straight, central, elastic, smooth
            %
            % v_p:  velocity of the puck
            % v_e:  velocity of the obstacle 
            
            % Divide velocity into tangential and normal velocity
            N = x_II - x_I;
            N = N / norm(N);
            T = [N(1); -N(2)];
            v0_I_T = (T' * v0_I) * T;
            v0_I_N = (N' * v0_I) * N;
            v0_II_T = (T' * v0_II) * T; 
            v0_II_N = (N' * v0_II) * N;
            
            % Calculate change in normal velocity
            v1_II_N = (2*v0_I_N*m_I - v0_II_N*m_I + v0_II_N*m_II)/(m_I + m_II);
            v1_I_N = (v0_I_N*m_I + v0_II_N*m_II - v1_II_N*m_II)/m_I; 
            
            % Calculate actual body velocities
            v1_I = v0_I_T + v1_I_N;
            v1_II = v0_II_T + v1_II_N;
            
            % Set positions
            diff = (r_I + r_II) - norm(x_II - x_I);
            if diff >= 0
                % TODO: This is a bit messy, maybe there is a better
                % solution
                x_I = x_I - N*diff/2;
                x_II = x_II + N*diff/2;
            end
        end
  	end
end
