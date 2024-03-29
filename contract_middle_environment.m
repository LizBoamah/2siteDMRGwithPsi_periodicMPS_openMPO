function middle_env = contract_middle_environment(M,H,M_dagger,numTensors)
% Initialize middle_env to an empty array
        
        % Compute middle environment
        middle_env = tensorprod(M{2}, H{2}, 2,2,"NumDimensionsA",3);
        middle_env = tensorprod(middle_env, M_dagger{2}, 4,2,"NumDimensionsA",5);

        % Contract remaining tensors to build middle environment
        for i = 3:(numTensors-1)
            temp_env = tensorprod(M{i}, H{i}, 2,2,"NumDimensionsA",3);
            temp_env = tensorprod(temp_env, M_dagger{i}, 4,2,"NumDimensionsA",5);
            middle_env = tensorprod(middle_env, temp_env, [2, 4, 6], [1, 3, 5]);
            middle_env = permute(middle_env, [1 4 2 5 3, 6]);
        end
    end