classdef Michalewicz
% The additive Michalewicz benchmark function
%
% Date: 15.08.2019
% Author: Michael Werner   

    properties(Access=public)
        num_vars;
        vars;
        minimize = true;
    end
    
    properties(Access=private)
       d = 10;
       m = 0.5;
    end
    
    methods
        function obj = Michalewicz()
            num_vars = obj.d;
            vars = [];
            for i=1:num_vars
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [0,pi])];
            end
            obj.vars = vars;
        end
        
        function s = call(obj, varargin)
            %CALL calculates the value of the Michalewicz function
            if length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            i = 1:obj.d;
            s = - sum(sin(varargin).* sin(i.*varargin.^2/pi)^(2*obj.m));
        end
    end
end