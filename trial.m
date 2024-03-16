% Construct the Hamiltonian
N = 6; % Number of lattice sites
t = 1; % Hopping parameter
Ho = construct_Hamiltonian(t, N);

% Compute the eigenvalues of the Hamiltonian
eigenValues = eig(Ho);

% Compute the expected eigenvalues from the cosine function
k = 0:N-1;
expectedEigenValues = -2 * t * cos(2 * pi * k / N);

% Sort the eigenvalues for proper visualization
eigenValues = sort(eigenValues);
expectedEigenValues = sort(expectedEigenValues);

% Plot the eigenvalues
figure;
hold on;
plot(eigenValues, 'o', 'DisplayName', 'Numerical Eigenvalues');
plot(expectedEigenValues, 'x', 'DisplayName', 'Expected Eigenvalues');
xlabel('State Index');
ylabel('Energy');
legend;
hold off;





































% % Define Parameters
% U = 1; % On-site interaction
% t = 1; % Hopping parameter
% N = 6; % Number of lattice sites
% 
% % Basis generation, 2^N states for N sites
% dim = 2^N;
% basis = de2bi((0:dim-1).');
% 
% % Initialization of the Hamiltonian
% H = zeros(dim, dim);
% 
% % Construct the Hamiltonian
% for i = 1:dim
%     state_i = basis(i, :);
%     for j = 1:dim
%         state_j = basis(j, :);
%         
%         % Hopping term - loop over each site except the last one due to OBC
%         for site = 1:(N - 1)
%             if state_i(site) ~= state_j(site) && state_i(site + 1) == state_j(site + 1)
%                 H(i, j) = H(i, j) - t;
%             end
%         end
%         
%         % Interaction term
%         for site = 1:(N - 1)
%             H(i, i) = H(i, i) + U * state_i(site) * state_i(site + 1);
%         end
%     end
% end
% 
% % Diagonalize the Hamiltonian
% [eigvec, eigval] = eig(H);
% 
% % Display Eigenvalues
% disp('Eigenvalues:');
% disp(diag(eigval));




% % Define Parameters
% U = 1; % On-site interaction
% t = 1; % Hopping parameter
% N = 6; % Number of lattice sites
% 
% % Basis generation, 2^N states for N sites
% dim = 2^N;
% basis = zeros(dim, N); % Initialize basis
% 
% for i = 0:(dim - 1)
%     binaryStr = dec2bin(i, N); % Convert to binary string with N digits
%     basis(i + 1, :) = str2double(regexp(num2str(binaryStr),'\d','match')); % Convert string to array of doubles
% end
% 
% % Initialization of the Hamiltonian
% H = zeros(dim, dim);
% 
% % Construct the Hamiltonian
% for i = 1:dim
%     state_i = basis(i, :);
%     for j = 1:dim
%         state_j = basis(j, :);
%         
%         % Hopping term - loop over each site except the last one due to OBC
%         for site = 1:(N - 1)
%             if state_i(site) ~= state_j(site) && state_i(site + 1) == state_j(site + 1)
%                 H(i, j) = H(i, j) - t;
%             end
%         end
%         
%         % Interaction term
%         for site = 1:(N - 1)
%             H(i, i) = H(i, i) + U * state_i(site) * state_i(site + 1);
%         end
%     end
% end
% 
% % Diagonalize the Hamiltonian
% [eigvec, eigval] = eig(H);
% 
% % Display Eigenvalues
% disp('Eigenvalues:');
% disp(diag(eigval));

