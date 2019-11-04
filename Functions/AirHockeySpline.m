classdef AirHockeySpline
% Date: 14.10.2019
% Author: Michael Werner   

    properties(Access=public)
        num_vars
        vars
        minimize = true;
    end
    
    properties(Access=private)
       InitialState
       SimulationParameter
       objectiveFunc
    end
    
    methods
        function obj = AirHockeySpline()
            mode = 1;
            switch mode
                case 1
                    obj.InitialState.x_puck = [0;0];
                    obj.InitialState.x_mallet = [0;-0.5];
                case 2
                    obj.InitialState.x_puck = [0;0];
                    obj.InitialState.x_mallet = [0.25;-0.75];
                case 3
                    obj.InitialState.x_puck = [-0.25;0];
                    obj.InitialState.x_mallet = [0;-0.5];
                case 4
                    obj.InitialState.x_puck = [-0.25;0];
                    obj.InitialState.x_mallet = [0.25;-0.75];
            end
            
            obj.InitialState.v_puck = [0;0];
            obj.InitialState.v_mallet = [0;0.0];
            obj.InitialState.goal_count = [0; 0];
            
            obj.SimulationParameter.dt = 0.01;
            obj.SimulationParameter.T = 5;
            
            taskAiming = AirHockeySim.TaskAimingSpline(obj.InitialState,obj.SimulationParameter);
            [vars, objectiveFunc] = taskAiming.CreateOptimizationTask();
            obj.vars = vars;
            obj.objectiveFunc = objectiveFunc;
            obj.num_vars = length(obj.vars);
        end
        
        function s = call(obj, varargin)
            s = feval(obj.objectiveFunc, varargin);
        end
        
        function s = visualize(obj, varargin)
           taskAiming = AirHockeySim.TaskAimingSpline(obj.InitialState,obj.SimulationParameter);
           s = taskAiming.evaluatePolicy(varargin, 1);
        end
    end
end


