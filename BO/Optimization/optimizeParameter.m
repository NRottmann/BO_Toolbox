function params = optimizeParameter(func, num_param)
    % Wrapper of fminunc that allows the user to group the parameters
    % num_param is a list containing the number of parameters per group

    % Get total number of hyperparameters     
    nvar = sum(num_param);

    % Use fminunc
    param_0 = 0.1*ones(nvar,1);
    options = optimoptions(@fminunc,'Display','off');
    param = fminunc(@callback,param_0,options);
    
    function f = callback(param)
       % split parameter into groups and store in cell array
       f = feval(func, groupParameters(param, num_param));
    end
        
    params = groupParameters(param, num_param);
end

function p = groupParameters(param, num_param)
   param_idx = 0;
   p = cell(length(num_param), 1);
   for group_idx=1:length(num_param)
      n = num_param(group_idx);
      p(group_idx) = {param(param_idx+1:param_idx+n)};
      param_idx = param_idx + n;
   end
end