classdef Rosenbrock
% The additive Rosenbrock benchmark function
% https://www.sfu.ca/~ssurjano/rosen.html
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
        function obj = Rosenbrock(d)
            obj.d = d;
            obj.num_vars = d;
            vars = [];
            for i=1:d
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [-5,10])];
            end
            obj.vars = vars;
        end
        
        function s = call(obj, varargin)
            %CALL calculates the value of the Rosenbrock function
            if  ~iscell(varargin) && length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            x = cell2mat(varargin);
            if isstruct(x)
               x = cell2mat(struct2cell(x));
            end
            x1 = x(1:end-1);
            x2 = x(2:end);
            s = sum(100*(x2-x1.^2).^2 + (1-x1).^2);
        end
    end
end