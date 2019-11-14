this folder contains the implemtation of typical benchmark functions and an air hockey simulation.

## How to add new functions
Each function is defined as a class object. It is required that the class has the public attribures vars and minimize. Additionally a public method call is required, to compute the function.
One can use the following template for new functions.
```matlab
classdef Func
% An example functions, calculating the sum of the parameters
    properties(Access=public)
        num_vars;
        vars;
        minimize = true;
    end
   
    methods
        function obj = Func(d)
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
            %CALL calculates the value of the function
            if length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            x = cell2mat(varargin);
            
            s = sum(x);
        end
    end
end
```
