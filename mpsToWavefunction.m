function psi = mpsToWavefunction(mps)
    % Convert an MPS to a conventional wavefunction (Psi)

    % Initialize the wavefunction with the first tensor of the MPS
    psi = reshape(mps{1}, [], size(mps{1}, 3)); % Reshape to a matrix

    % Iterate over the rest of the tensors in the MPS
    for i = 2:length(mps)
        tensor = mps{i};
        % Reshape tensor for matrix multiplication
        tensor = reshape(tensor, size(tensor, 1), []);

        % Contract tensor with wavefunction and reshape
        psi = psi * tensor;

        % Reshape psi to prepare for next contraction
        if i < length(mps)
            nextTensorDim = size(mps{i+1}, 1);
            psi = reshape(psi, [], nextTensorDim);
        end

    end

    % Reshape the final psi into a vector and normalize
    psi = reshape(psi, [], 1);
   
end


































% function Psi = mpsToWavefunction(M)
%     % Convert an MPS to a full wave function Psi.
%     % Inputs:
%     % M - Cell array representing the MPS, where each cell contains a tensor for a site.
%     %
%     % Output:
%     % Psi - Vector representing the full wave function.
% 
%     N = length(M); % Number of sites in the MPS
%     Psi = 1; % Initialize Psi as a scalar (will be expanded iteratively)
% 
%     for site = 1:N
%         tensor = M{site};
%         [rows, cols, height] = size(tensor);
%         Psi = reshape(Psi, [numel(Psi), 1, 1]); % Reshape Psi for tensor product
%         % Psi = ncon({Psi, tensor}, {[-(1:numel(Psi)), 1], [1, -numel(Psi)-1, -numel(Psi)-2]});
% %         Psi= tensorprod(Psi, tensor,2,1);
% %         Psi = reshape(Psi, [rows * cols * height, 1]); % Reshape back to vector
% %     end
% 
%     % Normalize Psi if needed
%     Psi = Psi / norm(Psi);
% end


% function psi = mpsToWavefunction(mps)
%     % Function to convert an MPS to a conventional wavefunction (Psi)
% 
%     % Initialize the wavefunction with the first tensor of the MPS
%     psi = reshape(mps{1}, [], 1); % Reshape the first tensor into a vector
% 
%     % Iterate over the rest of the tensors in the MPS
%     for i = 2:length(mps)
%         % Get the current tensor and reshape it into a matrix for the tensor product
%         tensor = reshape(mps{i}, size(mps{i}, 1), []);
% 
%         % Tensor product and contraction
%         psi = tensorprod(psi, tensor, 2, 1);
% 
%         % Reshape the resulting tensor back into a vector
%         psi = reshape(psi, [], 1);
%     end
% 
%     % Normalize Psi if necessary
%     psi = psi / norm(psi);
% end
























