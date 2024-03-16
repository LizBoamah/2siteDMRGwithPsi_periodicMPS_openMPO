function [isNormalized, isCanonical] = testMPS(M, dir, s)
    % Test if an MPS is normalized and in the correct canonical form
    % Inputs:
    % M - The MPS as a cell array of tensors
    % dir - The direction of canonicalization ('left', 'right', or 'mixed')
    % s - Orthogonality center (relevant for 'mixed' case)

    N = length(M); % Number of sites in the MPS

    % Check normalization
    psi = mpsToWavefunction(M); % Convert MPS to full wave function
    isNormalized = abs(norm(psi) - 1) < 1e-10; % Check if norm is approximately 1

    % Initialize canonical check
    isCanonical = true;

    switch dir
        case 'left'
            for l = 1:N-1
                isCanonical = isCanonical && checkLeftCanonical(M{l});
            end

        case 'right'
            for l = N:-1:2
                isCanonical = isCanonical && checkRightCanonical(M{l});
            end

        case 'mixed'
            for l = 1:s-1
                isCanonical = isCanonical && checkLeftCanonical(M{l});
            end
            for l = N:-1:s+2
                isCanonical = isCanonical && checkRightCanonical(M{l});
            end
    end
end

% Helper functions for left and right canonical checks
function isLeftCanonical = checkLeftCanonical(tensor)
    contracted = tensorprod(conj(tensor), tensor, [1, 3], [1, 3]);
    identityCheck = eye(size(contracted, 2));
    isLeftCanonical = all(abs(contracted(:) - identityCheck(:)) < 1e-10);
end

function isRightCanonical = checkRightCanonical(tensor)
    contracted = tensorprod(conj(tensor), tensor, [2, 3], [2, 3]);
    identityCheck = eye(size(contracted, 1));
    isRightCanonical = all(abs(contracted(:) - identityCheck(:)) < 1e-10);
end
