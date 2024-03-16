function isRightCanonical = checkRightCanonical(tensor)
    contracted = tensorprod(conj(tensor), tensor, [2, 3], [2, 3]);
    identityCheck = eye(size(contracted, 1));
    isRightCanonical = all(abs(contracted(:) - identityCheck(:)) < 1e-10);
end