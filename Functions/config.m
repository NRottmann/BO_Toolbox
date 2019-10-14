% Config file for the air hockey simulation environment
function out = config(caseString)
out = [];
switch caseString
    case 'AirHockeyInitialSetting'
        % States
        out.x_puck = [0; 0];
        out.x_mallet = [0; -0.5];
        out.v_puck = [0; 0];
        out.v_mallet = [0; 0];
        % Time step
        out.dt = 0.01;
        % Goal count
        out.goal_count = [0; 0];  
end
end