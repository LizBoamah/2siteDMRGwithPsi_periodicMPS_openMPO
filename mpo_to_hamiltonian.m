% function H = mpo_to_hamiltonian(MPO)
% % Converts MPO to Hamiltonian matrix
% % Input:
% %   MPO: cell array containing MPO tensors
% % Output:
% %   H: Hamiltonian matrix
% 
% N = length(MPO);
% d = size(MPO{1}, 2);
% 
% % Calculate the total Hilbert space dimension
% hilbert_dim = d^N;
% 
% % Initialize Hamiltonian matrix
% H = zeros(hilbert_dim);
% 
% % Loop over all basis states
% for state_dec = 0:(hilbert_dim-1)
%     state_bin = dec2bin(state_dec, N) - '0';
%     bra = state_bin + 1;
% 
%     for new_state_dec = 0:(hilbert_dim-1)
%         new_state_bin = dec2bin(new_state_dec, N) - '0';
%         ket = new_state_bin + 1;
% 
%         coef = 1;
%         bond_idx1 = ones(1, N);
%         bond_idx2 = ones(1, N);
%         for i = 1:N
%             D1 = size(MPO{i}, 1);
%             D2 = size(MPO{i}, 4);
%             local_coef = 0;
%             for b1 = 1:D1
%                 for b2 = 1:D2
%                     local_coef = local_coef + MPO{i}(b1, bra(i), ket(i), b2) * bond_idx1(i) * bond_idx2(i);
%                 end
%             end
%             coef = coef * local_coef;
%             bond_idx1(i) = mod(bond_idx1(i), D1) + 1;
%             bond_idx2(i) = mod(bond_idx2(i), D2) + 1;
%         end
% 
%         H(state_dec+1, new_state_dec+1) = coef;
%     end
% end
% end
function H_mpo = mpo_to_hamiltonian(MPO)

    N = numel(MPO); % assuming MPO is a cell array
    if N < 2
        error('MPO must contain at least two elements.');
    end
    nextDim = 6; % Start with 6 as mentioned 6-D when i=3
    temp = tensorprod(MPO{1}, MPO{2}, 4, 1);
    for i = 3:N
        temp1 = tensorprod(temp, MPO{i}, nextDim, 1);
        temp = temp1; % Update temp to hold the result of the current iteration
        nextDim = nextDim + 2; % Increase dimensionality by 2 for the next iteration
    end   
    % Get the size of the original tensor
       sz = size(temp1);

    % Number of dimensions of the tensor
      numDims = length(sz);

    % Generate permutation order dynamically
       permOrder = [1:2:numDims, 2:2:numDims];
    

    % Permute the tensor
    B = permute(temp1, permOrder);

    
    % Calculate the products of the sizes of odd and even dimensions
    s1 = prod(sz(1:2:numDims));
    s2 = prod(sz(2:2:numDims));

    % Reshape the permuted tensor
    H_mpo = reshape(B, [s1, s2]);

end



































% function H = mpo_to_hamiltonian(MPO)
% % Converts MPO to Hamiltonian matrix
% % Input:
% %   MPO: cell array containing MPO tensors
% % Output:
% %   H: Hamiltonian matrix
% 
% N = length(MPO);
% d = size(MPO{1}, 2);
% 
% % Calculate the total Hilbert space dimension
% hilbert_dim = d^N;
% 
% % Initialize Hamiltonian matrix
% H = zeros(hilbert_dim);
% 
% % Loop over all basis states
% for state_dec = 0:(hilbert_dim-1)
%     state_bin = dec2bin(state_dec, N) - '0';
%     bra = state_bin + 1;
%     
%     for new_state_dec = 0:(hilbert_dim-1)
%         new_state_bin = dec2bin(new_state_dec, N) - '0';
%         ket = new_state_bin + 1;
%         
%         coef = 1;
%         for i = 1:N
%             D1 = size(MPO{i}, 1);
%             D2 = size(MPO{i}, 4);
%             for bond_idx1 = 1:D1
%                 for bond_idx2 = 1:D2
%                     coef = coef * MPO{i}(bond_idx1, bra(i), ket(i), bond_idx2);
%                 end
%             end
%         end
%         
%         H(state_dec+1, new_state_dec+1) = coef;
%     end
% end
% end



















% function H = mpo_to_hamiltonian(MPO)
% % Converts MPO to Hamiltonian matrix
% % Input:
% %   MPO: cell array containing MPO tensors
% % Output:
% %   H: Hamiltonian matrix
% 
% N = length(MPO);
% d = size(MPO{1}, 2);
% 
% % Calculate the total Hilbert space dimension
% hilbert_dim = d^N;
% 
% % Initialize Hamiltonian matrix
% H = zeros(hilbert_dim);
% 
% % Loop over all basis states
% for state_dec = 0:(hilbert_dim-1)
%     state_bin = dec2bin(state_dec, N) - '0';
%     bra = state_bin + 1;
%     
%     for new_state_dec = 0:(hilbert_dim-1)
%         new_state_bin = dec2bin(new_state_dec, N) - '0';
%         ket = new_state_bin + 1;
%         
%         coef = 1;
%         for i = 1:N
%             coef = coef * MPO{i}(bra(i), state_bin(i)+1, ket(i), new_state_bin(i)+1);
%         end
%         
%         H(state_dec+1, new_state_dec+1) = coef;
%     end
% end
% end









