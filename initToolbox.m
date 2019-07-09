function initToolbox()
% Function to add allfolders to the matlab path
paths = genpath ( pwd );
exclude = {'.','..','.git'};

% Get all the paths between the PATHSEP character
pat = ['([^', pathsep, ']+)', pathsep, '?'];
paths = cellfun(@(x) x{1}, regexp(paths, pat, 'tokens'), 'UniformOutput', false);

% Filter out unwanted directories
mask = 0;
for i=1:length(exclude)
ex = cell2mat(exclude(1, i));
mask = mask + contains(paths, ex);
end
paths = paths(~mask);

% Finally add the paths
addpath(paths{:});
end