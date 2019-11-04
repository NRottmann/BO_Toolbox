classdef AirHockey
% Date: 29.08.2019
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
        function obj = AirHockey()
            obj.InitialState.x_puck = [0;0];
            obj.InitialState.v_puck = [0;0];
            obj.InitialState.x_mallet = [0;-0.5];
            obj.InitialState.v_mallet = [0;0.0];
            obj.InitialState.goal_count = [0; 0];
            
            obj.SimulationParameter.dt = 0.01;
            obj.SimulationParameter.T = 5;
            
            taskAiming = AirHockeySim.TaskAiming(obj.InitialState,obj.SimulationParameter);
            [vars, objectiveFunc] = taskAiming.CreateOptimizationTask();
            obj.vars = vars;
            obj.objectiveFunc = objectiveFunc;
            obj.num_vars = length(obj.vars);
        end
        
        function s = call(obj, varargin)
            s = feval(obj.objectiveFunc, varargin);
        end
        
        function s = visualize(obj, varargin)
           s = obj.taskAiming.evaluatePolicy(varargin, 1);
        end
    end
end


