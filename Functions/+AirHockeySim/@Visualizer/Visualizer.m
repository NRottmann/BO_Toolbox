% The MIT License (MIT)
% Copyright (c) 2017 Elmar Rueckert
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

% Visualizer for the Air Hockey Simulator
% 
% Author: Nils Rottmann
% E-Mail: Rottmann@rob.uni-luebeck.de
%
% Date: 19.03.2019

classdef Visualizer < handle
    properties (Access = public)
         window;                 % The overall window used for visualization
    end
    properties (Access = private)
        classKineticModel;      % The model which shall be used for visualization
        
        w_win = 700;            % Width and height of the visualization window
        h_win;
        
        r_puck;                 % Game data for visualization
        r_mallet;
        r_center
        w_field;
        l_field;
        w_goal;
        l_goal;
        
        x_field;                % Rectangular objects
        y_field;
        x_edge;
        y_edge;
        x_goal_t;
        y_goal_t;
        x_goal_b;
        y_goal_b;
        
        x_puck;                 % Circular objects
        y_puck;
        x_center;
        y_center;
        x_mallet;
        y_mallet;

        draw_phi = 0:pi/50:2*pi;    % defines how many points are used for drawing the circles
        
        table;                  % The table struct which holds all visualizationd data
    end
    methods (Access = public)
        % CONSTRUCTOR
        function obj = Visualizer(kineticModel)
            % Allocate variables
            obj.classKineticModel = kineticModel;
            
            % Get game data from the kinetic model
            data = kineticModel.getGameData();
            obj.r_puck = data.r_puck;
            obj.r_mallet = data.r_mallet;
            obj.r_center = data.r_puck * 2;
            obj.w_field = data.w_field;         
            obj.l_field = data.l_field;
            obj.w_goal = data.w_goal;
            obj.l_goal = data.w_goal * 0.6;
            
            % Calculating the position of the rectangular objects
            obj.x_field = obj.w_field/2 * [-1 -1 1 1 -1];
            obj.y_field = obj.l_field/2 * [-1 1 1 -1 -1];
            obj.x_edge = (obj.w_field+0.12)/2 * [-1 -1 1 1 -1];     % the edges are 0.06 m wide
            obj.y_edge = (obj.l_field+0.12)/2 * [-1 1 1 -1 -1];
            obj.x_goal_t = (obj.w_goal)/2 * [-1 -1.2 1.2 1 -1];     % position for top goal area
            obj.y_goal_t = obj.l_field/2 + 0.06 * [0 1 1 0 0];
            obj.x_goal_b = (obj.w_goal)/2 * [-1 -1.2 1.2 1 -1];     % position for bottom goal area
            obj.y_goal_b = -obj.l_field/2 - 0.06 * [0 1 1 0 0];
            
            % Calculating the position of the circular objects
            data = kineticModel.getPositionData();
            obj.x_puck = data.x_puck(1) + obj.r_puck * cos(obj.draw_phi);
            obj.y_puck = data.x_puck(2) + obj.r_puck * sin(obj.draw_phi);
            obj.x_center = obj.r_center * cos(obj.draw_phi);
            obj.y_center = obj.r_center * sin(obj.draw_phi);
            obj.x_mallet = data.x_mallet(1) + obj.r_mallet * cos(obj.draw_phi);
            obj.y_mallet = data.x_mallet(2) + obj.r_mallet * sin(obj.draw_phi);
            
            % Define window size
            obj.h_win = obj.w_win * obj.l_field / obj.w_field;
        end
        
        % Initialize the visualization window with the table
        function obj = InitializeVisualization(obj)
            % Calculating the position of the circular objects
            data = obj.classKineticModel.getPositionData();
            obj.x_puck = data.x_puck(1) + obj.r_puck * cos(obj.draw_phi);
            obj.y_puck = data.x_puck(2) + obj.r_puck * sin(obj.draw_phi);
            obj.x_center = obj.r_center * cos(obj.draw_phi);
            obj.y_center = obj.r_center * sin(obj.draw_phi);
            obj.x_mallet = data.x_mallet(1) + obj.r_mallet * cos(obj.draw_phi);
            obj.y_mallet = data.x_mallet(2) + obj.r_mallet * sin(obj.draw_phi);
            
            % Drawing the main window
            obj.window = figure('unit','pixels','numbertitle','off','menubar','none','resize','off',...
                    'position',[0 0 obj.w_win obj.h_win],'name','Air Hockey Simulator');
            movegui('center');
            % Drawing the "container" which contains all the gameplay objects 
            obj.table.container = axes('parent',obj.window,'unit','pixels');
            set(obj.table.container,'dataaspectratio',[1 1 1]);
            hold on
            % Drawing the edge around the hockey table
            obj.table.edge = fill(obj.x_edge,obj.y_edge,'b','Facecolor',[0 0 1],'parent',obj.table.container);
            hold on
            % Drawing the playing area
            obj.table.field = fill(obj.x_field,obj.y_field,'b','facecolor', [1 1 1],'edgecolor','none','parent',obj.table.container);
            hold on
            % Drawing the center circle
            obj.table.centercircle = plot(obj.x_center, obj.y_center, 'b','color',[0 0 1],'linewidth',4,'parent',obj.table.container);
            hold on  
            % Drawing the bottom goal crease
            obj.table.creaseb = plot(obj.w_goal*[-1/2 -1/2 1/2 1/2], ...
                        [-obj.l_field/2 -obj.l_field/2 + obj.l_goal -obj.l_field/2 + obj.l_goal -obj.l_field/2],...
                           'color',[0 0 1],'linewidth',2,'parent',obj.table.container);
            hold on
            % Drawing the top goal crease
            obj.table.creaset = plot(obj.w_goal*[-1/2 -1/2 1/2 1/2], ...
                        [obj.l_field/2 obj.l_field/2 - obj.l_goal obj.l_field/2 - obj.l_goal obj.l_field/2],...
                           'color',[0 0 1],'linewidth',2,'parent',obj.table.container);
            % Drawing the centerline
            obj.table.centerline = plot([-obj.w_field/2 obj.w_field/2],[0 0],'color',[0 0 1],'linewidth',4,'parent',obj.table.container);
            hold on
            % Drawing the puck
            obj.table.puck = fill(obj.x_puck,obj.y_puck,'r','parent',obj.table.container);
            hold on
            % Drawing the user's mallet
            obj.table.mallet = fill(obj.x_mallet, obj.y_mallet,'k','parent',obj.table.container);
            hold on
            obj.table.dot = fill(0 + obj.r_mallet/5 * cos(obj.draw_phi), 0 + obj.r_mallet/5 * sin(obj.draw_phi),'g','parent',obj.table.container);
            hold on
            obj.table.dotv = fill(0 + obj.r_mallet/5 * cos(obj.draw_phi), 0 + obj.r_mallet/5 * sin(obj.draw_phi),'g','parent',obj.table.container);
            hold on
            obj.table.t = fill(0 + obj.r_mallet/5 * cos(obj.draw_phi), 0 + obj.r_mallet/5 * sin(obj.draw_phi),'g','parent',obj.table.container);
            hold on
            % Drawing the top goal area
            obj.table.goal_t = fill(obj.x_goal_t,obj.y_goal_t,'b','Facecolor',[0 0 0],'edgecolor', 'none','parent',obj.table.container);
            hold on
            % Drawing the bottom goal area
            obj.table.goal_b = fill(obj.x_goal_b,obj.y_goal_b,'b','Facecolor',[0 0 0],'edgecolor', 'none','parent',obj.table.container);
            % Setting the display limits and hiding the axes
            set(obj.table.container,'xlim',[-obj.w_field/2 - 0.05, obj.w_field/2 + 0.05]);
            set(obj.table.container,'ylim',[-obj.l_field/2 - 0.05, obj.l_field/2 + 0.05]);
            axis off;
            drawnow;
        end
        
        % Update the position of the puck and the mallet
        function obj = UpdateVisualization(obj)
            % Get position of the puck and mallet
            data = obj.classKineticModel.getPositionData();
            obj.x_puck = data.x_puck(1) + obj.r_puck * cos(obj.draw_phi);
            obj.y_puck = data.x_puck(2) + obj.r_puck * sin(obj.draw_phi);
            obj.x_mallet = data.x_mallet(1) + obj.r_mallet * cos(obj.draw_phi);
            obj.y_mallet = data.x_mallet(2) + obj.r_mallet * sin(obj.draw_phi);
            % Drawing the puck
            delete(obj.table.puck);
            obj.table.puck = fill(obj.x_puck,obj.y_puck,'r','parent',obj.table.container);
            hold on
            % Drawing the user's mallet
            delete(obj.table.mallet);
            obj.table.mallet = fill(obj.x_mallet, obj.y_mallet,'k','parent',obj.table.container);
            drawnow
        end
        
        % Update the position of the puck and the mallet
        function obj = UpdateVisualization_advanced(obj, x1, y1, vx1, vy1, q)
            s = sqrt(vx1 * vx1 + vy1 + vy1);
            vx1 = vx1/s * 0.25;
            vy1 = vy1/s * 0.25;
            % Get position of the puck and mallet
            data = obj.classKineticModel.getPositionData();
            obj.x_puck = data.x_puck(1) + obj.r_puck * cos(obj.draw_phi);
            obj.y_puck = data.x_puck(2) + obj.r_puck * sin(obj.draw_phi);
            obj.x_mallet = data.x_mallet(1) + obj.r_mallet * cos(obj.draw_phi);
            obj.y_mallet = data.x_mallet(2) + obj.r_mallet * sin(obj.draw_phi);
            % Drawing the puck
            delete(obj.table.puck);
            obj.table.puck = fill(obj.x_puck,obj.y_puck,'r','parent',obj.table.container);
            hold on
            % Drawing the user's mallet
            delete(obj.table.mallet);
            obj.table.mallet = fill(obj.x_mallet, obj.y_mallet,'k','parent',obj.table.container);
            hold on
            delete(obj.table.dot);
            obj.table.dot = fill(x1 + obj.r_mallet/5 * cos(obj.draw_phi), y1 + obj.r_mallet/5 * sin(obj.draw_phi),'g','parent',obj.table.container);
            delete(obj.table.dotv);
            obj.table.dotv = fill(x1 + obj.r_mallet/7 * cos(obj.draw_phi) + vx1, y1 + obj.r_mallet/7 * sin(obj.draw_phi) + vy1,'g','parent',obj.table.container);
            hold on
            delete(obj.table.t)
            obj.table.t = plot(q(1,:), q(2,:), 'b','color',[0 0 1],'parent',obj.table.container);
            drawnow
        end
    end
end