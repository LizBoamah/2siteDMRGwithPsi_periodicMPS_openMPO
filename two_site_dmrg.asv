function [full_sweep_energy, new_energy_values, M, E_exact,occupation_number_old,occupation_number_new] = two_site_dmrg(N, bd, U,max_dims, D,D_periodic,mu, t, NrEl, max_sweeps)

% TWO_SITE_DMRG - Perform two-site DMRG for the Hubbard model
%
% N - number of sites
% bd - bond dimension
% D - maximum bond dimension
% U - onsite interaction strength
% t - hopping parameter
% max_sweeps - maximum number of sweeps
% tol - convergence tolerance
% energy - ground state energy
% M - MPS representation (cell array of tensors)

 


 % Initialize a basis states to form a initial wave function
   B=mps_form_full_basis_new(bd,N,[0 1]',NrEl);

 % Calculate the dimension vector for each site
 dimvec = repelem(bd,N);

 % Calculate the total dimension of the full Hilbert space
 Dim = prod(dimvec);

 % Calculate the dimension of the subspace
  DimSubspace = size(B,1);
  % fill with random values
     C=rand(DimSubspace,1);
 
     % fill Psi
        Psi=zeros(Dim,1);
     ind = zeros(DimSubspace,1);
     %   % Loop over each state in the subspace
     for i=1:DimSubspace
         % Calculate the index in the full Hilbert space
          Bi=num2cell(B(i,:));
         % Bi=fliplr(num2cell(B(i,:)));
         ind(i)=sub2ind(dimvec,Bi{:});
         % Assign the amplitude from C to the corresponding position in Psi
            Psi(ind(i)) = C(i);
     end%for i

    % Normalize Psi
       Psi = Psi / norm(Psi);
       

    
% Hubbard Hmailtonian (2nd quantatised-full space)  for comparsion with
% Hmatrix
 Ho = construct_Hamiltonian(t,U, N);
 [~,~] = exact_diagonalization(Ho);

% Construct the MPO for the Hubbard Hamiltonian(hard coded)
 H = hubbard_mpo_site(t,U,mu, N);
 %convert the MPO to hamiltonian matrix
   Hmatrix = mpo_to_hamiltonian(H);
  % Diagonalise to find the lowest exact energy
   [~,E_exact] = exact_diagonalization(Hmatrix);
% % Set direction: right, left, mixed
    dir = 'mixed';

 % Initialize the MPS representation
    M = mps_canonical(Psi,bd, N, dir,1,max_dims);
      M= normalizeMPS(M, dir,1);
     % calcuate  the conjugate of each MPS tensor in the cell array
       numTensors = length(M);  %  M is a 1D cell array
        M_dagger = cell(1, numTensors);
            for i = 1:numTensors
                M_dagger{i} = conj(M{i});
            end
               M_dagger= normalizeMPS(M_dagger, dir,1);
            
        % Some test functions to check correct of MPS after conversion
           % Psi_old = mpsToWavefunction(M);
         [~, ~] = testMPS(M, dir, 1); % checks normalization and  if MPS is canonicalised.
        numberMPO = occupation_mpo_site(N);
        occupation_number_old = occupation_number_calculation(M,numberMPO, M_dagger);% Occupation number calculation before dmrg
    
        % all_symmetric = is_MPS_symmetric(M); % checks the symmetry of the MPS.
      


% DMRG step begins here . 
% Start the sweeping process
for sweep = 1
    full_sweep_energy =  zeros(1, 2*sweep);
    E_left = zeros(sweep,1);
    % Sweep from left to right
    for s = 1:N-1 
        % Environment is contracted for easy computations
       [left_env,middle_env, right_env] = contract_environments(M, M_dagger, H, s, N);

        % Combine the environments and the two-site tensor to build the effective Hamiltonian
        effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env,s,N,middle_env);
    
        % Diagonalize the effective Hamiltonian
         [eig_vec, eig_val] = eig((effective_H+effective_H')/2);
         energy_new1 = sort(diag(eig_val),'ascend');
         energy_new1 = energy_new1(1);
         idx = find(diag(eig_val) == eig_val, 1);
          psi = eig_vec(:, idx);
           E_left(sweep,s) = energy_new1;
          full_sweep_energy(s,2*sweep-1) = energy_new1;
           
      % Reshape psi_matrix and perform SVD
        psi_matrix = reshape(normalize(psi), [size(left_env,1)*size(H{s},2) , size(H{s+1}, 2)*size(right_env,1) ]);
        [U, S, V] = svd(psi_matrix, 'econ');
            % % Truncate the matrices
            % if size(S, 1) > max_dims
            %     U = U(:, 1:max_dims);
            %     S = S(1:max_dims, 1:max_dims);
            %     V = V(:, 1:max_dims);
            % end

          % % Update M{s} and M{s+1}
          % Update M{s} and M{s+1}
          M{s} = reshape(U, [size(left_env,1), size(H{s}, 2), size(S, 1)]);
          M{s+1} = reshape(S*V', [size(S, 2), size(H{s+1}, 2), size(right_env, 1)]);

        % Update the conjugate MPS tensor M_
            M_dagger{s} = conj(M{s});
            M_dagger{s+1} = conj(M{s+1});
  
   end

occupation_number_new = occupation_number_calculation(M,numberMPO, M_dagger); % occupation number calculation after dmrg


% Two-site update for the (N,1) edge after the first sweep
 % Special handling for connecting the first and last sites
    s = N; % Setting s = N to handle the connection between site N and 1
    [left_env, middle_env, right_env] = contract_environments(M, M_dagger, H, s, N);
    effective_H = build_effective_hamiltonian(left_env, H{s}, H{1}, right_env,s,N, middle_env); % Notice H{1} for the last site

    % Diagonalize the effective Hamiltonian
    [eig_vec, eig_val] = eig((effective_H+effective_H')/2);
    [energy, ground_state_index] = min(diag(eig_val));
    psi = eig_vec(:, ground_state_index);

    % Store the energy
     % E_left(s,) = energy;
     full_sweep_energy(s,2*sweep-1) = energy;

    % Reshape psi and perform SVD
    psi_matrix = reshape(normalize(psi), [size(middle_env, 2) * size(M{s}, 2), size(M{1}, 2) * size(middle_env, 1)]);
    [U, S, V] = svd(psi_matrix, 'econ');

    % Truncate to maximum bond dimension D
    if size(S, 1) > D_periodic
        U = U(:, 1:D_periodic);
        S = S(1:D_periodic, 1:D_periodic);
        V = V(:, 1:D_periodic);
    end
   % Update M and M_dagger for the current edge
     % Update M and M_dagger for the (N,1) edge
    M{s} = reshape(U, [ size(middle_env, 2), size(M{s}, 2),size(S, 1)]);
    M{1} = reshape(S*V',[size(S, 2), size(M{1}, 2) ,size(middle_env, 1)]);

    % Update the conjugate MPS tensor M_dagger
    M_dagger{s} = conj(M{s});
    M_dagger{1} = conj(M{1});

end

        % Continue with subsequent sweeps (from 2 to max_sweeps) with periodic MPS
        % Initialise new energy values for storing the energies with
        % periodic boundary conditions
         % new_energy_values= zeros(max_sweeps,N);
         new_energy_values = zeros(N, 2*max_sweeps);
          new_E_left =zeros(max_sweeps,N);
          new_E_right = zeros(max_sweeps,N);

    for sweep = 1:max_sweeps
        % forward and backward sweep code including the (N,1) edge
        % Forward sweep
        for s = 1:N
        % Contract the environments
        if s < N
            [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, s, N);
            % Compute the new effective Hamiltonian
            updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{s}, H{s+1}, updated_right_env, updated_middle_env);
        else
            % Special case for the edge between the last and the first site
            [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, s, N);
            updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{s}, H{1}, updated_right_env, updated_middle_env);
        end

                 % Diagonalize the effective Hamiltonian
                 [eig_vec, eig_val] = eig((updated_effective_H + updated_effective_H') / 2);
                 energy_new1 = sort(diag(eig_val),'ascend');
                 energy_new1 = energy_new1(1);
                 idx = find(diag(eig_val) == eig_val, 1);
                  psi = eig_vec(:, idx);
                  new_E_left(sweep,s) = energy_new1;
                  new_energy_values(s,2*sweep) = energy_new1;



                % Reshape psi and perform SVD
                if s<N
                     psi_matrix = reshape(psi, [size(M{s}, 1)*size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
                else
                     % Reshape psi and perform SVD for edge(1)=N
                     psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{1}, 2) * size(M{1}, 3)]);
                end

                [U, S, V] = svd(psi_matrix, 'econ');
                % coefficients = diag(S);

            % Truncate to maximum bond dimension D
                if size(U, 2) > D
                    U = U(:, 1:D);
                    S = S(1:D, 1:D);
                    V = V(:, 1:D);
                end
                % Update M and M_dagger using the ground state wavefunction
                % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
                  if s == 1
                    % Update M{edge} and M{edge+1}
                     M{s} = reshape(U, [size(updated_right_env,2), size(M{s}, 2), size(S, 1)]);
                     M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(updated_right_env, 1)]);
                elseif s == N-1
                    % Update M{edge} and M{edge+1}
                    M{s} = reshape(U, [size(updated_left_env, 2), size(M{s}, 2), size(S, 1)]);
                    M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(updated_left_env,1)]);
                  elseif s ==  N
                     M{s} = reshape(U, [size(updated_middle_env, 2), size(M{s}, 2), size(S, 1)]);
                     M{1} = reshape(S*V', [size(S, 2), size(M{1}, 2), size(updated_middle_env, 1)]);
                  else
                       M{s} = reshape(U, [size(updated_left_env, 2), size(M{s}, 2),size(S, 1)]);
                       M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2),size(updated_right_env, 1)]);
                 end

                % Update the conjugate MPS tensor M_dagger
                if s<N
                    M_dagger{s} = conj(M{s});
                    M_dagger{s+1} = conj(M{s+1});
                else
                    M_dagger{s} = conj(M{s});
                    M_dagger{1} = conj(M{1});
                end
        end


    end
% Flatten E_left and E_right for easier processing
new_E_left_flat = reshape(new_E_left.', 1, []);
new_E_right_flat = reshape(new_E_right.', 1, []);

% Calculate the error
new_error_left = abs(new_E_left_flat - E_exact);
new_error_right = abs(new_E_right_flat - E_exact);

num_sweeps = size(new_E_left, 1);
num_sites = size(new_E_left, 2);

% Prepare microstep indices for plotting
microstep_indices = 1:(num_sweeps * num_sites);

% Plotting the errors
figure;
plot(microstep_indices, new_error_left, '-o', 'DisplayName', 'Error Left Sweep');
hold on;
plot(microstep_indices, new_error_right, '-x', 'DisplayName', 'Error Right Sweep');
xlabel('Microstep');
ylabel('Error (Absolute Difference from Exact Energy)');
title('DMRG Convergence Analysis');
legend;
grid on;























%         disp(lowest_energy);

% % Check convergence after a full sweep
%     if sweep > 1  % Skip the first sweep as we don't have two energy values yet
%         energy_diff = abs(energy_values(end) - energy_values(end-1));
%         if energy_diff < tol
%             break;  % Stop the sweeping process if the solution has converged
%         end
%     end
% end
% 
% if converged
%     fprintf('DMRG converged after %d sweeps.\n', sweep);
% else
%     fprintf('DMRG did not converge within the maximum number of sweeps.\n');
% end






% for s = N-1:-1:1
%        [left_env,middle_env, right_env] = contract_environments(M, M_dagger, H, s, N);
% 
%         % Combine the environments and the two-site tensor to build the effective Hamiltonian
%         effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env,s,N,middle_env);
% 
%         % Diagonalize the effective Hamiltonian
%          [eig_vec, eig_val] =  eig((effective_H+effective_H')/2);
%          energy_new2 = sort(diag(eig_val),'ascend');
%          energy_new2 = energy_new2(1);
%          idx = find(diag(eig_val) == eig_val, 1);
%           psi = eig_vec(:, idx);
%           E_right(sweep,s) =  energy_new2;
%           full_sweep_energy(s,2*sweep) = energy_new2;
% 
% 
%         psi_matrix = reshape(normalize(psi), [size(left_env,3)*size(H{s},2) , size(H{s+1}, 2)*size(right_env,1) ]);
%         [U, S, V] = svd(psi_matrix, 'econ');
%         %   % Truncate the matrices
%         % if size(S, 1) > max_dims
%         %     U = U(:, 1:max_dims);
%         %     S = S(1:max_dims, 1:max_dims);
%         %     V = V(:, 1:max_dims);
%         % end
% 
%          % Update M{s} and M{s+1}
%           M{s} = reshape(U*S, [size(left_env,1), size(H{s}, 2), size(S, 1)]);
%           M{s+1} = reshape(V', [size(S, 2), size(H{s+1}, 2), size(right_env, 1)]);
% 
%         % Update the conjugate MPS tensor M_
%         M_dagger{s} = conj(M{s});
%         M_dagger{s+1} = conj(M{s+1});   
















%         % Backward sweep
%         for s = N:-1:1
%              if s < N
%                 [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, s, N);
%                 % Compute the new effective Hamiltonian
%                 updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{s}, H{s+1}, updated_right_env, updated_middle_env);
%             else
%                 % Special case for the edge between the last and the first site
%                 [updated_left_env, updated_middle_env, updated_right_env] = updated_contract_environments(M, M_dagger, H, s, N);
%                 updated_effective_H = build_updated_effective_hamiltonian(updated_left_env, H{s}, H{1}, updated_right_env, updated_middle_env);
%             end
% 
%                 % Diagonalize the effective Hamiltonian
%                  [eig_vec, eig_val] =eig((updated_effective_H + updated_effective_H') / 2);
% 
%                  energy_new2 = sort(diag(eig_val),'ascend');
%                  energy_new2 = energy_new2(1);
%                  idx = find(diag(eig_val) == eig_val, 1);
%                   psi = eig_vec(:, idx);
%                   new_E_right(sweep,s) = energy_new2;
%                   new_energy_values(N+1-s,2*sweep-1) = energy_new2;
%                   % energy_values(N+1-s,2*sweep-1) = energy_new2
% 
%                   % Reshape psi and perform SVD
%                 if s<N
%                      psi_matrix = reshape(psi, [size(M{s}, 1)*size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%                 else
%                      % Reshape psi and perform SVD for edge(1)=N
%                      psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{1}, 2) * size(M{1}, 3)]);
%                 end
% 
%             [U, S, V] = svd(psi_matrix, 'econ');
%             % Truncate to maximum bond dimension D
%             if size(U, 2) > D
%                 U = U(:, 1:D);
%                 S = S(1:D, 1:D);
%                 V = V(:, 1:D);
%             end
% 
%               % Update M and M_dagger using the ground state wavefunction
%                 % This typically involves reshaping psi and possibly performing an SVD to maintain the MPS form
%                   if s == 1
%                     % Update M{s} and M{s+1}
%                      M{s} = reshape(U*S, [size(updated_right_env,2), size(M{s}, 2), size(S, 1)]);
%                      M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(updated_right_env, 1)]);
%                 elseif s == N-1
%                     % Update M{s} and M{s+1}
%                     M{s} = reshape(U*S, [size(updated_left_env, 2), size(M{s}, 2), size(S, 1)]);
%                     M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(updated_left_env,1)]);
%                   elseif s ==  N
%                      M{s} = reshape(U*S, [size(updated_middle_env, 2), size(M{s}, 2), size(S, 1)]);
%                      M{1} = reshape(V', [size(S, 2), size(M{1}, 2), size(updated_middle_env, 1)]);
%                   else
%                       M{s} = reshape(U*S, [size(updated_left_env, 2), size(M{s}, 2),size(S, 1)]);
%                       M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2),size(updated_right_env, 1)]);
%                  end
% 
%            % Update the conjugate MPS tensor M_dagger
%                 if s<N
%                     M_dagger{s} = conj(M{s});
%                     M_dagger{s+1} = conj(M{s+1});
%                 else
%                     M_dagger{s} = conj(M{s});
%                     M_dagger{1} = conj(M{1});
%                 end
% 
%         end



















% for s = N-1:-1:1
% %      % Initialize the MPS representation
% %       M = mps_canonical(Psi,bd, N, dir,s);
% %     % Calculate the complex conjugate of each MPS tensor in the cell array
% %     M_ = cellfun(@conj, M, 'UniformOutput', false);
% 
%     % Perform the necessary updates and contractions for right-to-left sweep
%     if s == N-1
%         % Contract environment from the left
%         left_env = 1;
%         for j = 1:s-1
%             left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%         end
%         right_env = 1;
%     elseif s == 1
%         % Contract environment from the right
%         right_env = 1;
%         for j = N:-1:s+2
%             right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%         end
%         left_env = 1;
%     else
%         % Contract environment from the left
%         left_env = 1;
%         for j = 1:s-1
%             left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%         end
%         
%         % Contract environment from the right
%         right_env = 1;
%         for j = N:-1:s+2
%             right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%         end
%     end
% 
%     % Combine the environments and the two-site tensor to build the effective Hamiltonian
%     effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
% 
%     % Diagonalize the effective Hamiltonian
%     [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
%     energy_new = (min(diag(eig_val)));
%     idx = find(diag(eig_val) == eig_val, 1);
%     psi = eig_vec(:, idx);
%     E_right = [E_right; energy_new];
%    
% 
% %   fprintf('Sweep: %d, Sites: (%d, %d), Energy: %f\n', sweep, s, s+1, energy_new);
% %         
%     if s == N-1
%         % Reshape psi and perform SVD for s=N
%         psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     elseif s == 1
%         % Reshape psi and perform SVD for s=1
%         psi_matrix = reshape(psi, [size(M{s}, 2), size(M{s+1}, 1) * size(M{s+1}, 3)]);
%     else
%         % Reshape psi and perform SVD for other values of s
%         psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     end
% 
%     [U, S, V] = svd(psi_matrix, 'econ');
% 
%    if s == 1
%         % Update M{s} and M{s+1}
%         M{s} = reshape(U * S, [1, size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%     elseif s == N-1
%         % Update M{s-1} and M{s}
%         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(S * V', [size(S, 2), size(M{s+1}, 2), 1]);
%     else
%         % Update M{s} and M{s+1}
%         M{s} = reshape(U * S, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%     end
% 
%     % Update the conjugate MPS tensor M_
%     M_{s} = conj(M{s});
%     M_{s+1} = conj(M{s+1});
% end
%   % Store energy value after completing one full sweep
%         E_full_sweep = (E_left(end) + E_right(end));
%         energy_values = [energy_values; E_full_sweep];
% 
% 
% end






































%         if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [1, size(M{s}, 2), size(S, 1)]);
%             M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%     %         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     %         M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%              M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%         end

%  % Perform the necessary updates and contractions for left-to-right sweep
%         if s == 1
%             % Contract environment from the right
%             right_env = 1;
%             for j = N:-1:s+2
%                 right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%             end
%             left_env = 1;
%         elseif s == N-1
%             % Contract environment from the left
%             left_env = 1;
%             for j = 1:s-1
%                 left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%             end
%             right_env = 1;
%         else
%             % Contract environment from the left
%             left_env = 1;
%             for j = 1:s-1
%                 left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%             end
%             
%             % Contract environment from the right
%             right_env = 1;
%             for j = N:-1:s+2
%                 right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%             end
%        end
% 
% 

% if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             end
% 
% 
%         % Update the conjugate MPS tensor M_
%         M_{s} = conj(M{s});
%         M_{s+1} = conj(M{s+1});




%         if s == 1
%             % Initialize the MPS representation
%             M = mps_canonical(Psi,bd, N, dir, s);
%             % Calculate the complex conjugate of each MPS tensor in the cell array
%             M_ = cellfun(@conj, M, 'UniformOutput', false);
%         
%         end
%        [left_env, right_env] = contract_environments(M, M_, H, s, N);
% 
%         % Combine the environments and the two-site tensor to build the effective Hamiltonian
%         effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
%     
%         % Diagonalize the effective Hamiltonian
%         [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
%         energy_new = (min(diag(eig_val)));
%         idx = find(diag(eig_val) == eig_val, 1);
%         psi = eig_vec(:, idx);
%           
%          E_left = [E_left; energy_new];
%          
% 
%             % Calculate the new energy after completing a left-to-right sweep
%        
%       % Reshape psi and perform SVD for s=1
%          if s == 1
%             psi_matrix = reshape(psi, [size(M{s}, 2)* size(M{s},1), size(M{s+1}, 1) * size(M{s+1}, 2)]);
%         elseif s == N-1
%             % Reshape psi and perform SVD for s=N
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         else
%             % Reshape psi and perform SVD for other values of s
%             psi_matrix = reshape(psi, [size(M{s}, 2) * size(M{s}, 3), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         end
%     
%           [U, S, V] = svd(psi_matrix, 'econ');
% 
%        if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [ size(M{s}, 2), size(M{s}, 1), size(S, 1)]);
%             M{s+1} = reshape(S*V', [ size(M{s+1}, 3),size(M{s+1}, 1), size(M{s+1}, 2) ]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%              M{s} = reshape(U, [size(M{s}, 2), size(M{s}, 3), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 2)]);
%        end
% 
%         % Update the conjugate MPS tensor M_
%         M_{s} = conj(M{s});
%         M_{s+1} = conj(M{s+1});
% 
%        
% 
%       
%   
% 
% 
%      end
%     








%          if s == 1
%             psi_matrix = reshape(psi, [size(M{s}, 2)*size(M{s},3), size(M{s+1}, 1) * size(M{s+1}, 2)]);
%         elseif s == N-1
%             % Reshape psi and perform SVD for s=N
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         else
%             % Reshape psi and perform SVD for other values of s
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         end
%     



%        if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [size(S, 1),size(M{s}, 2), size(M{s}, 3)]);
%             M{s+1} = reshape(S*V', [size(M{s+1}, 1), size(M{s+1}, 2),size(S, 2)]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%     %         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     %         M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%              M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%        end




%         % Truncate singular values
%             truncation_idx = min(size(S, 1), D);
%             S = S(1:truncation_idx, 1:truncation_idx);
%             U = U(:, 1:truncation_idx);
%             V = V(:, 1:truncation_idx);
% 
%           if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%           end













% function[energy_new, M] = two_site_dmrg(N, bd, D, U, s, t, max_sweeps, tol)
% % TWO_SITE_DMRG - Perform two-site DMRG for the Hubbard model
% % 
% % N - number of sites
% % bd - bond dimension
% % D - maximum bond dimension
% % U - onsite interaction strength
% % t - hopping parameter
% % max_sweeps - maximum number of sweeps
% % tol - convergence tolerance
% 
% %energy - ground state energy
% %M - MPS representation (cell array of tensors)
% 
% %Construct the Hubbard Hamiltonian 
% % Ho = construct_Hamiltonian(t, U, N);
% 
% %Initialize a MPS 
% mps = init_random_mps(N, bd, D);
% 
% %Set direction: right, left, mixed
% dir = 'mixed';
% 
% %Construct the MPO for the Hubbard Hamiltonian
% % H = MPOcompress(Ho, bd, N);
%  H= hubbard_mpo_site(U, t, N, bd, D);
% %Initialize the MPS representation
% M = mps_canonicalM(mps, N,dir, s);
% 
% %Calculate the complex conjugate of each MPS tensor in the cell array
% 
%  M_ = cellfun(@conj, M,'UniformOutput', false);
% %Initialize variables for sweeping
%   energy_old = inf;
% 
% %Start the sweeping process
% for sweep = 1:max_sweeps
%     %Contract environment from the left
%     left_env = 1;
%     for j = 1:s-1
%         left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%     end
% 
%     %Contract environment from the right
%     right_env = 1;
%     for j = N:-1:s+2
%         right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%     end
% 
%     %Combine the environments and the two-site tensor to build the effective Hamiltonian
%     effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env);
% 
%    % Diagonalize the effective Hamiltonian
%     [eig_vec, eig_val] = eig(effective_H);
%     energy_new = real(min(diag(eig_val)));
%     idx = find(diag(eig_val) == energy_new, 1);
%     psi = eig_vec(:, idx);
% 
%    % Reshape psi and perform SVD
%     psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     [U, S, V] = svd(psi_matrix, 'econ');
% 
% %     %Truncate the bond dimension
% %     if size(S, 2) > D
% %         S = S(:, 1:D);
% %         U = U(:, 1:D);
% %         V = V(1:D, :);
% %     end
%     % Print the energy at the current pair of sites
%     fprintf('Sweep: %d, Sites: (%d, %d), Energy: %f\n', sweep, s, s+1, energy_new);
% 
% 
%     %Update MPS tensors
%     M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     M{s+1} = reshape(S*V, [size(S, 1), size(M{s+1}, 2), size(M{s+1}, 3)]);
% 
%  % Check for convergence
%     energy_diff = abs(energy_new - energy_old);
%     if energy_diff < tol
%         
%     end
%     energy_old = energy_new;
% end
% 
%    
% end
% 
% 



% % Check convergence
%     energy_diff = abs(energy_new - energy_old);
%     if energy_diff < tol
% %         converged = true;
%         break;
%     else
%         energy_old = energy_new;
%     end
% end
% energies = [energies, energy];
% 



%         % Truncate singular values
%             truncation_idx = min(size(S, 1), D);
%             S = S(1:truncation_idx, 1:truncation_idx);
%             U = U(:, 1:truncation_idx);
%             V = V(:, 1:truncation_idx);
% 
%           if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%           end

































