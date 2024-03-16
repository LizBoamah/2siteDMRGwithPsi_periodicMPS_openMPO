function MPO = hubbard_mpo_site( t,U, mu,N)
% Constructs the MPO tensor for the Hubbard Hamiltonian in the second
% quantized form for a one-dimensional chain of N sites.
% Input:
%   U: onsite interaction strength
%   t: hopping parameter
%   N: number of sites
%   d: bond dimension
% Output:
%   MPO: cell array containing MPO tensors

% Create MPO for Hubbard Hamiltonian
MPO = cell(1, N);
% d = 2;
% D= 10;

% Identity and fermion creation/annihilation matrices
I = eye(2);
c = [0, 1; 0, 0];
c_dag = c';
n = c_dag * c;
ph = diag([1,-1]);

% Construct MPO tensors
for i = 1:N
    if i == 1
        W = zeros(1, 2, 2, 5);
        W(1, :, :, 1) = I;
        W(1, :, :, 2) = -t *c_dag* ph;
        W(1, :, :, 3) = t * c*ph;
        W(1, :, :, 4) = U*n;
        W(1, :, :, 5) = mu *n;
    elseif i == N
        W = zeros(5, 2, 2, 1);
        W(1, :, :, 1) = mu*n;
        W(2, :, :, 1) = c;
        W(3, :, :, 1) = c_dag;
        W(4, :, :, 1) = n;
        W(5, :, :, 1) = I;
       
    else
        W = zeros(5, 2, 2, 5);
        W(1, :, :, 1) = I;
        W(1, :, :, 2) = -t *c_dag *ph;
        W(1, :, :, 3) = t *c *ph;
        W(1, :, :, 4) = U *n;
        W(1, :, :, 5) = mu*n;
        W(2, :, :, 5) = c;
        W(3, :, :, 5) = c_dag;
        W(4, :, :, 5) = n;
        W(5, :, :, 5) = I;

    end
    MPO{i} = W;
end
end


























