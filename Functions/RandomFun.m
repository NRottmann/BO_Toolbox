classdef RandomFun
    % TODO
    
    properties
        X;
        Y;
        Z;
        x_f;
        z_f;
        X_idx_F;
    end
    
    methods
        function obj = RandomFun()
            load('randomData.mat','X','Y','Z','x_f','z_f','X_idx_F');
            obj.X = X;
            obj.Y = Y;
            obj.Z = Z;
            obj.x_f = x_f;
            obj.z_f = z_f;
            obj.X_idx_F = X_idx_F;
        end
        
        function zq = paramFun(obj,xq,yq)
            zq = interp2(obj.X,obj.Y,obj.Z,xq,yq);
        end
        
        function xf = featureFun(obj,xq,yq)
            % Compute Euclidean distances:
            distances = sqrt(sum(bsxfun(@minus, obj.X_idx_F, [xq, yq]).^2,2));
            % Find the smallest distance and use that as an index into B:
            [~,idx] = min(distances);
            xf = obj.x_f(idx);
        end
    end
end

