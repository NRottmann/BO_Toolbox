classdef NeuralNet
    %NEURALNET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=public)
        num_param
    end
    
    properties(Access=private)
        architecture
    end
    
    methods(Access = private)
       function n_weights = getnumweights(obj)
        % Computes the number of weights for a given architecture.
            n_weights = 0;
            for layer_idx=1:length(obj.architecture)-1
                n_weights = n_weights + (obj.architecture(layer_idx) + 1) * obj.architecture(layer_idx +1);
            end
       end 
       
       function f = probagate(obj, x, weights)
           %PROBAGATE probagtes the input x through the neural net.
           if(size(x, 1) ~= obj.architecture(1))
             error("Input size does not match input layer size")
           end

           y = x';
           w = weights;
           for layer_idx=1:length(obj.architecture)-1
             W = w(1:obj.architecture(layer_idx) * obj.architecture(layer_idx+1));
             W = reshape(W, obj.architecture(layer_idx), obj.architecture(layer_idx+1));
             b = w(obj.architecture(layer_idx) * obj.architecture(layer_idx+1) + 1 : ...
                 (obj.architecture(layer_idx) + 1) * obj.architecture(layer_idx+1));
             a =  y * W + b;
             y =  1.0 ./ (1 + exp(-a));
             
             w = w((obj.architecture(layer_idx) + 1) * obj.architecture(layer_idx+1):end); 
           end
           f = y';
       end
    end
    
    methods(Access=public)
        function obj = NeuralNet(D, d, inputVars, varargin)
            %NEURALNET Construct an instance of this class
            %   Detailed explanation goes here
            
            % Default values
            defaultargs = {'num_hidden', 20}; 
            params = setargs(defaultargs, varargin);
            obj.architecture = [D, params.num_hidden, d];
            obj.num_param = obj.getnumweights();
        end
               
        function f = getfeature(obj, x, param)
           if length(param) ~= obj.num_param
               error(strcat("The given number of parameters does not match the requested size of ", num2str(obj.num_param), "."));
           end
           f = obj.probagate(x, param');
        end
               
        function bounds = getbounds(obj)
            %GETBOUNDS returns the bounds of the feature space           
            % activation function 1/(1+exp(-a)) returns values from 0 to 1
            bounds = [zeros(obj.architecture(end), 1), ones(obj.architecture(end), 1)];
        end
    end
end