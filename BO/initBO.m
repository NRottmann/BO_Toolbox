function initBO()
% Function to add allfolders to the matlab path

% Get the folder names
nameFolds = cell(1,1);
nameFolds{1} = 'AcqFunctions';

% Get the path of your directory
full = mfilename('fullpath');
name = mfilename;
base = full(1:end-length(name));

% Construct the code which should be inserted into the startup script
for i=1:1:length(nameFolds)
    str_add = strcat([base nameFolds{i}]);
    addpath(str_add);
end
end