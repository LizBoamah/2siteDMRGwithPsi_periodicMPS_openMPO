function H = MPOcompress(Ho, d, N)
    ml = 1;

    % Initialize MPO - tensor at site n
    H = cell(1, N);

    % Loop from first to N-1 site
    for l = 1:N - 1
        % Reshape Ho to prepare for SVD
        W = reshape(Ho, [ml * d * d, numel(Ho) / (ml * d * d)]);

        % Perform SVD
        [U, S, V] = svd(W, 'econ');

        % Determine new bond dimension after SVD
        new_ml = size(U, 2);

        % Reshape U to make it a rank-4 tensor of shape [ml, d, d, new_ml]
        U = reshape(U, [ml, d, d, new_ml]);

        % Prepare next Ho for following iteration
        Ho = S * V';
        
        % Store the tensor
        H{l} = U;

        % Update the bond dimension for next iteration
        ml = new_ml;
    end

    % Handle the last site
    H{N} = reshape(Ho, [new_ml, d, d, 1]);
end

% 
% function H = MPOcompress(Ho,d, N)
% 
%     ml = 1;
% 
%     % Initialize MPO - tensor at site n
%     H = cell(1, N);
%     %--- beginning to N-1
%     for l = 1:N - 1
%         W = reshape(Ho, [ml * d^d, d^(2*(N - l))]); 
%         [U, S, V] = svd(W, 'econ');
%         new_ml = size(U, 2);
% 
%         U = reshape(U, [ml, d, d, new_ml]);
%         V = S * V';
% 
%         H{l} = U;
% 
%         Ho = V;
%         ml = new_ml;
%     end
% 
%     %---last mode
%     H{N} = reshape(Ho, [new_ml, d, d, 1]);
%     
% 
%     
% end
