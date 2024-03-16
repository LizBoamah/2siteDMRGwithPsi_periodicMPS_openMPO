function normalizedMPS = normalizeMPS(M, dir,s)
    % Function to normalize a Matrix Product State (MPS) in canonical form

    numTensors = length(M);  % Number of tensors in the MPS
    normalizedMPS = M;  % Copy the MPS

    if strcmp(dir, 'right')
        % Normalize the first tensor for right-canonical form
        firstTensor = M{1};
        normFactor = sqrt(sum(conj(firstTensor(:)) .* firstTensor(:)));
        normalizedMPS{1} = firstTensor / normFactor;

    elseif strcmp(dir, 'left')
        % Normalize the last tensor for left-canonical form
        lastTensor = M{numTensors};
        normFactor = sqrt(sum(conj(lastTensor(:)) .* lastTensor(:)));
        normalizedMPS{numTensors} = lastTensor / normFactor;

    elseif strcmp(dir, 'mixed')
        % Find and normalize the center tensor for mixed canonical form
        % Assuming the center tensor is known or can be identified
        centerTensor = M{s};  % Define centerTensorIndex based on your MPS structure
        normFactor = sqrt(sum(conj(centerTensor(:)) .* centerTensor(:)));
        normalizedMPS{s} = centerTensor / normFactor;
    end
end
