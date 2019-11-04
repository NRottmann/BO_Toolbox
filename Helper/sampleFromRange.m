function samples = sampleFromRange(dim, num_samples, data)
%SAMPLEFROMRANGE Generate uniformly distributed sample distribution
if ~isa(data(1), 'double')
    % convert list of 'optimizableVariable' into range matrix
    range = reshape([data(:).Range], 2, [] )';
else
    range = data;
end
samples = rand(dim, num_samples) ... % random matrix, uniform in [0,1]
    .* repmat(range(:, 2) - range(:,1), 1, num_samples)... % scale each dimension
    + range(:, 1); % translate in to corrent range
end

