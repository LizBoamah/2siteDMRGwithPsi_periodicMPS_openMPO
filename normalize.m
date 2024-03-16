function x_normalized = normalize(x)
    norm_x = norm(x); % Compute the 2-norm of x
    x_normalized = x / norm_x; % Divide by its norm to normalize
end

