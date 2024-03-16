function testfor_two_site_dmrg

% [lowest_energy,energy_values, M,E_exact] = two_site_dmrg(N, bd, U, t, max_sweeps, tol)
%  [lowest_energy,energy_values, M,E_exact] = two_site_dmrg(6, 2, 1, 1, 2, 20);


N=5;
d=2;
NrEl=2;
D=7; % bond dimension


dimvec = repelem(d,N);
Dim = prod(dimvec);

B=mps_form_full_basis_new(d,N,[0 1]',NrEl);
DimSubspace = size(B,1);
C=rand(DimSubspace,1);


%% fill Psi
Psi=zeros(Dim,1);
for i=1:DimSubspace
    Bi=num2cell(B(i,:));
    ind=sub2ind(dimvec,Bi{:});
    Psi(ind) = C(i);
end%for i

Psi = Psi/norm(Psi);

%% check the number of electrons
id=eye(d);
% n1=[0,0;0,1];
n_tot_op = zeros(Dim);
n_expval=zeros(1,N);
for l=1:N
%     ni=1;
%     for ll=1:N
%        if ll==l
%            ni=kron(ni,n1);
%        else
%            ni=kron(ni,id);
%        end%if
%    end%for
%    n_tot_op = n_tot_op + ni;
    ni = annihilate_op(l,N);
    ni = ni' * ni;

    n_expval(l) = Psi' * ni * Psi;

    n_tot_op = n_tot_op + ni;
end%for

sum(n_expval)
Psi' * n_tot_op * Psi % this should be equal to NrEl

%% Psi to Mps with truncation
M = mps_canonical_left_trunc(Psi, d, N, D);

%% check left-canonical
% M = mps_canonical(Psi,d, N, 'left',1);
for l=1:N
    Ml=M{l};
    % Sum=zeros(size(Ml,3));
    % Mlq = squeeze(Ml(:,q,:));
    Sum = tensorprod(Ml,conj(Ml),[1,2],[1,2]);
    max( abs(Sum-eye(size(Ml,3))),[], "all" )
end%for

%% check Psi and Mps values

ii = randi(DimSubspace,1);
bb = num2cell(B(ii,:)); 
ind=sub2ind(dimvec,Bi{:});
Psi(ind)
C_Mps = squeeze( M{1}(:,bb{1},:) )';
for l=2:N
    C_Mps = C_Mps * squeeze( M{l}(:,bb{l},:) );
end%for
C_Mps

% for i=1:DimSubspace;
%     Bi=num2cell(B(i,:));
%     ind=sub2ind(dimvec,Bi{:});
%     C_Psi = Psi(ind)
%     % C_M   = 1)
%     tensorprod()
% end%for i

% Psi_full = M{1};
% for l=2:N
%     Psi_full = tensorprod(Psi_full,M{l},[1+l],[1])
% end%for

end%function testfor_two_site_dmrg


function M = mps_canonical_left_trunc(Psi, d, N, D)
% this shoul be a separate function

% mr=size(Psi,1);
% case 'left'
M = cell(1,N);
ml = 1;
for l = 1:N - 1
    W = reshape(Psi, [ml * d, d ^ (N - l)]);
    [U, S, V] = svd(W, 'econ');
    new_ml = size(U, 2); new_ml = min(new_ml,D);

    U_trunc = U(:,1:new_ml);
    S_trunc = S(1:new_ml,1:new_ml);
    V_trunc = V(:,1:new_ml);

    M{l} = reshape(U_trunc, [ml, d, new_ml]);
    Psi = S_trunc * V_trunc';
    ml = new_ml;
end
%---last mode
M{N} = reshape(Psi, [ml, d, 1]);

end%function mps_canonical_trunc



% %% check right-canonical
% M = mps_canonical(Psi,d, N, 'right',1);
% for l=1:N
%     Ml=M{l};
%     S=zeros(size(Ml,1));
%     for q=1:d
%         Mlq = squeeze(Ml(:,q,:));
%         S=S+Mlq*Mlq';
%     end%for q
%     S
% end%for
% 
% %% mixed-canonical
% M = mps_canonical(Psi,d, N, 'mixed',1);
% 
