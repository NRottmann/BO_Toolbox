classdef Schwefel
% The additive Schwefel benchmark function
% https://www.sfu.ca/~ssurjano/schwef.html
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
        function obj = Schwefel(d)
            obj.d = d;
            obj.num_vars = d;
            vars = [];
            for i=1:d
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [-500,500])];
            end
            obj.vars = vars;
        end
        
        function s = call(obj, varargin)
            %CALL calculates the value of the Schwefel function 
            
            if ~iscell(varargin) && length(varargin) ~= obj.d
               error(strcat('The number of input parameters (',...
                            num2str(length(varargin)), ...
                            ') does not match the number of dimension: ',...
                            num2str(obj.d)));
            end
            x = cell2mat(varargin);
            if isstruct(x)
               x = cell2mat(struct2cell(x));
            end
            s = 418.9829*obj.d - sum(x.*sin(sqrt(abs(x))));
        end
    end
end