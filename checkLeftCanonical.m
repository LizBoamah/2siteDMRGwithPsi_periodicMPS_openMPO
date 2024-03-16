function isLeftCanonical = checkLeftCanonical(tensor)
    contracted = tensorprod(conj(tensor), tensor, [1, 3], [1, 3]);
    identityCheck = eye(size(contracted, 2));
    isLeftCanonical = all(abs(contracted(:) - identityCheck(:)) < 1e-10);
end