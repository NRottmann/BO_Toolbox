classdef Ackley
% The Ackley benchmark function
% https://www.sfu.ca/~ssurjano/ackley.html
%
% Date: 15.08.2019
% Author: Michael Werner   

    properties(Access=public)
        num_vars
        vars
        minimize = true;
    end
    
    properties(Access=private)
       d = 10;
       a = 20;
       b = 0.2;
       c = 2*pi;
    end
    
    methods
        function obj = Ackley(d)
            obj.d = d;
            obj.num_vars = d;
            vars = [];
            for i=1:d
               vars = [vars, optimizableVariable(strcat('x', num2str(i)),...
                                                 [-5,5])];
            end
            obj.vars = vars;
        end
        
        function s = call(obj, varargin)
            %CALL calculates the value of the Ackley function          
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
            s = -obj.a * exp(-obj.b*sqrt(1/obj.d * sum(x.^2))) ...
                - exp(1/obj.d*sum(cos(obj.c*x))) ...
                + obj.a + exp(1);
        end
    end
end