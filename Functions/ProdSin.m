classdef ProdSin
% The non-additive Product of Sines benchmark function
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
        function obj = ProdSin()
            num_vars = obj.d;
            vars = [];
            for i=1:num_vars
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [0,2*pi])];
            end
            obj.vars = vars;
        end
        
        function p = call(obj, varargin)
            %CALL calculates the value of the Product of Sines function
            if length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            p = sin(varargin(1)) * prod(sin(varargin));
            
            if obj.minimize
                p = -p;
            end
        end
    end
end