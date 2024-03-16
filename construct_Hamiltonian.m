function Ho = construct_Hamiltonian(t, U, N)
    % Constructing a tridiagonal hopping matrix
    t_mat = diag(-t * ones(1, N-1), 1) + diag(-t * ones(1, N-1), -1);
    
    Ho = 0;
    for k = 1:N-1
        if U == 0
            Ho = Ho + (t_mat(k, k+1) * ((creation_op(k, N) * annihilate_op(k+1, N) + creation_op(k+1, N) * annihilate_op(k, N))));
        else
            Ho = Ho + t_mat(k, k+1) * (creation_op(k, N) * annihilate_op(k+1, N) + creation_op(k+1, N) * annihilate_op(k, N)) ...
                 + U * (num_op(k, N) * num_op(k+1, N));
                
          % Ho = Ho + t_mat(k, k+1)*(creation_op(N, N) * annihilate_op(1, N) +  creation_op(1, N)* annihilate_op(N, N) );
%         creation_opup(N, N) * annihilate_opup(1, N) +  creation_opup(1, N)* annihilate_opup(N, N));
    
        end

    end
     % % Making it a PBC(periodic boundary condition)
     %         Ho= Ho +(creation_op(N, N) * annihilate_op(1, N) +  creation_op(1, N)* annihilate_op(N, N) );
end

