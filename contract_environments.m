function [left_env,middle_env, right_env] = contract_environments(M, M_, H, s, N)
% Special case where edge is between the last and first sites
    if s == N 
         middle_env = contract_middle_environment(M,H,M_,N);
        left_env = 1;
        right_env = 1;

    elseif s == 1
        % Contract environment from the right
        right_env = 1;
        for j = N:-1:s+2
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
        left_env = 1;
        middle_env = 1;

    elseif s == N-1
        % Contract environment from the left
        left_env = 1;
        for j = 1:s-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end
        right_env = 1;
        middle_env = 1;
    else
        % Contract environment from the left
        left_env = 1;
        for j = 1:s-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end

        % Contract environment from the right
        right_env = 1;
        for j =  N:-1:s+2
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
        middle_env =1;
    end
end


    
