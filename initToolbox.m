function initToolbox()
% Function to add allfolders to the matlab path

% Get the folder names
d = dir;
isub = [d(:).isdir];                    % returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..','.git'})) = [];

% Get the path of your directory
full = mfilename('fullpath');
name = mfilename;
base = full(1:end-length(name));

% Construct the code which should be inserted into the startup script
for i=1:1:length(nameFolds)
    str_add = strcat([base nameFolds{i}]);
    addpath(str_add);
end

% Initialize subfolders
initGP

end