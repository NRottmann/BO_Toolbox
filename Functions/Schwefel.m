classdef Schwefel
% The additive Schwefel benchmark function
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
    end
    
    methods
        function obj = Schwefel()
            num_vars = obj.d;
            vars = [];
            for i=1:num_vars
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [-500,500])];
            end
            obj.vars = vars;
        end
        
        function s = call(obj, varargin)
            %CALL calculates the value of the Schwefel function
            if length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            s = 418.9829*obj.d - sum(varargin.*sin(sqrt(abs(varargin))));
        end
    end
end