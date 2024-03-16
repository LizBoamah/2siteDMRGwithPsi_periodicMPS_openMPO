function [M] = mps_canonical(Psi,d, N, dir,s, max_dim)

    ml = 1;
    m_mix =1;

    % Initialise MPS - tensor at site n
    M = cell(1,N);

    switch dir

        case 'left'
            for l = 1:N - 1
                W = reshape(Psi, [ml * d, d ^ (N - l)]); 
                [U, S, V] = svd(W, 'econ');
                % Truncate if necessary
                new_ml = min(size(U, 2), max_dim);
                U = U(:, 1:new_ml);
                S = S(1:new_ml, 1:new_ml);
                V = V(:, 1:new_ml);

                U = reshape(U, [ml, d, new_ml]);
                V = S * V';
                M{l} = U(:, :, :);
                Psi = V;
                ml = new_ml;
            end
            Psi = reshape(Psi, [ml, d, 1]);
            M{N} = Psi(:,:,:);

        case 'right'
            for l=N:-1:2
                W = reshape(Psi, d^(l-1), ml*d);
                [U, S, V] = svd(W, 'econ');
                % Truncate if necessary
                new_ml = min(size(V, 2), max_dim);
                U = U(:, 1:new_ml);
                S = S(1:new_ml, 1:new_ml);
                V = V(:, 1:new_ml);

                V = ctranspose(V);
                V = reshape(V, [new_ml, d, ml]);
                M{l} = V(:,:,:);
                Psi = U*S;
                ml = new_ml;
            end
            Psi = reshape(Psi, [1, d, ml]);
            M{1} = Psi(:,:, :);

        case 'mixed'
            % left canonical from 1 to s-1
            for l = 1:s-1
                W = reshape(Psi, [ml * d, d ^ (N - l)]);
                [U, S, V] = svd(W, 'econ');
                % Truncate if necessary
                new_ml = min(size(U, 2), max_dim);
                U = U(:, 1:new_ml);
                S = S(1:new_ml, 1:new_ml);
                V = V(:, 1:new_ml);

                M = reshape(U, [ml, d, new_ml]);
                V = S * V';
                M{l} = M(:, :, :);
                Psi = V;
                ml = new_ml;
            end

            % right canonical from N to s+1
            for l=N:-1:s+1
                W = reshape(Psi, d^(l-1), m_mix*d);
                [U, S, V] = svd(W, 'econ');
                % Truncate if necessary
                new_mlx = min(size(V, 2), max_dim);
                U = U(:, 1:new_mlx);
                S = S(1:new_mlx, 1:new_mlx);
                V = V(:, 1:new_mlx);

                V = ctranspose(V);
                V = reshape(V, [new_mlx, d, m_mix]);
                U = U*S;
                M{l} = V(:,:,:);
                Psi = U;
                m_mix = new_mlx;
            end

             %  othorgonal center Sites s
             W = reshape(Psi, [ml*d, m_mix]);
             [U, S, V] = svd(W, 'econ');
            % Truncate if necessary
            new_ml = min(size(S, 1), max_dim);
            U = U(:, 1:new_ml);
            S = S(1:new_ml, 1:new_ml);
            V = V(:, 1:new_ml);
            V= ctranspose(V);
            V= reshape(U*S, [ml, d, new_ml]);
            M{s} = V(:, :, :);






            % % Sites s and s+1
            % W = reshape(Psi, [ml*d, d*m_mix]);
            % [U, S, V] = svd(W, 'econ');
            % % Truncate if necessary
            % new_ml = min(size(U, 2), max_dim);
            % U = U(:, 1:new_ml);
            % S = S(1:new_ml, 1:new_ml);
            % V = V(:, 1:new_ml);
            % 
            % U = reshape(U, [ml, d, new_ml]);
            % M{s} = U(:, :, :);
            % V = reshape(S*V', [new_ml, d, m_mix]);
            % M{s+1} = V(:, :, :);

    end
end

% function [M] = mps_canonical(Psi,d, N, dir,s)
% 
%     ml = 1;
%     m_mix =1;
% 
%     % Initialise MPS - tensor at site n
%     M = cell(1,N);
% 
%   switch dir
% 
%         case 'left'
%             %--- beginning to N-1
%             for l = 1:N - 1
%                 W = reshape(Psi, [ml * d, d ^ (N - l)]); 
%                 [U, S, V] = svd(W, 'econ');
%                 new_ml = size(U, 2);
% 
%                 U = reshape(U, [ml, d, new_ml]);
%                 V = S * V';
%                M{l} = U(:, :, :);
%                 Psi = V;
%                 ml = new_ml;
%             end
%             %---last mode
%             Psi = reshape(Psi, [ml, d, 1]);
% %             for qq = 1:d
% %                 M{N}{qq} = permute(reshape(Psi(:, qq), [ml, 1, d]), [1, 3, 2]);
%                  M{N} = Psi(:,:,:);
% %             end
% 
% 
%          case 'right' % chk Ok
%             %--- N to 2
%             for l=N:-1:2
%               W=reshape(Psi, d^(l-1),ml*d);
%               [U,S,V]=svd(W, 'econ');
%               V = ctranspose(V);
%               new_ml = size(V, 1);
%               V = reshape(V, [new_ml, d, ml]);
%               M{l} = V(:,:,:);
%               Psi = U*S;
%               ml = new_ml;
%             end
%              %--- first mode
%             Psi = reshape(Psi, [1, d, ml]);
%               M{1} = Psi(:,:, :);
% 
%          case 'mixed' %chk ok
%              %left canonical from 1 to s-1
%                 for l = 1:s-1
%                     W = reshape(Psi, [ml * d, d ^ (N - l)]);
%                     [U, S, V] = svd(W, 'econ');
%                     new_ml = size(U, 2);
%                     M = reshape(U, [ml, d, new_ml]);
%                     V = S * V';
%                     M{l} = M(:, :, :);
%                     Psi = V;
%                     ml = new_ml;
%                 end
% 
% 
%                 % right canonical from N to s+1
%                 for l=N:-1:s+2
%                         W=reshape(Psi, d^(l-1),m_mix*d);
%                       [U,S,V]=svd(W, 'econ');
%                        V = ctranspose(V);
%                       new_mlx = size(V, 1);
%                       V = reshape(V, [new_mlx, d, m_mix]);
%                       U = U*S;
%                        M{l}= V(:,:,:);
%                        Psi = U;
%                       m_mix = new_mlx; 
% 
%                 end 
% 
%             % Sites s and s+1
%              W = reshape(Psi, [ml*d, d*m_mix]);
%              [U, S, V] = svd(W, 'econ');
%                new_ml = size(U, 2);
%                U = reshape(U, [ml, d, new_ml]);
%                M{s} = U(:, :, :);
%             V = reshape(S*V', [new_ml, d, m_mix]);
%             M{s+1} = V(:, :, :);
% 
% 
% 
% 
% 
%  end
% 
% 
%  % % Sites s and s+1
%  %            W = reshape(Psi, [ml * d, new_mlx]);
%  %            [U, S, V] = svd(W, 'econ');
%  % 
%  %               new_mlx = size(U, 1);
%  %               U = reshape(U*S, [ml, d, new_mlx]);
%  %              % U = reshape(U, [ml, new_mlx, new_mlx]);
%  %            M{s} = U(:, :, :);
%  %             new_mlx = size(V', 2);
%  %            V = reshape(V', [new_mlx, d, m_mix]);
%  %            M{s+1} = V(:, :, :);
%  % 
