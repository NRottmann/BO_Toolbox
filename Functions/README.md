this folder contains the implementation of typical benchmark functions and an air hockey simulation.

## How to use functions
Assumed we have a function called ```Func```. Then you can access the function by it attributes ```num_vars``` and ```vars```. Additionally you can compute the function value, using the ```call``` method of this class.
```matlab
f = Func(10);    % create 10 dimensional version of the function
v = f.vars;      % access data about the variables of the function (type is optimizableVariable)
n = f.num_vars;  % access the number of variables
s = f.call(1:10) % comput function value for input (1,...,10)
```

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
