function numberMPO = occupation_mpo_site(N)
% Constructs the MPO tensor for the Hubbard Hamiltonian in the second
% quantized form for a one-dimensional chain of N sites.
% Input:
%   U: onsite interaction strength
%   t: hopping parameter
%   N: number of sites
%   d: bond dimension
% Output:
%   MPO: cell array containing MPO tensors
% 
% Create MPO for Hubbard Hamiltonian
numberMPO = cell(1, N);
% d = 2;
% D= 10;

% Identity and fermion creation/annihilation matrices
% Define local operators in the spinless fermion basis: |0>, |1>
I = eye(2);
c = [0, 1; 0, 0]; % annihilation operator
c_dag = c';       % creation operator
n = c_dag * c;    % number operator
% ph = diag([1, -1]);
mu = 1;

% Construct MPO tensors
for i = 1:N
    if i == 1
        W = zeros(1, 2, 2, 2);
        W(1, :, :, 1) = I;
        W(1, :, :, 2) = mu*n;
    elseif i == N
        W = zeros(2,2,2,1);
        W(1, :, :, 1) = mu*n;
        W(2, :, :, 1) = I;
    else
        W = zeros(2, 2, 2, 2);
        W(1, :, :, 1) = I;
        W(1, :, :, 2) = mu*n;
        W(2, :, :, 2) = I;
    end
    numberMPO{i} = W;
 

end








