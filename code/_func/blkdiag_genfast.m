function Q = blkdiag_genfast(A, n)

% Performs fast n-times block-diagonalization of numeric (r x c) array A.
% Does the same job as MATLAB's function BLKDIAG, with the exception that
% the user does not need to manually specify repeated inputs of the array A
% n-times. This is useful when n is not known a priori -- e.g. when n
% varies in a loop.
%
% See also PERMUTE, RESHAPE, BLKDIAG

Q = zeros( [ size(A), n, n] );
Q(:, :, 1:n+1:n^2) = repmat(A, [1 1 n]);
Q = permute(Q, [1 3 2 4]);
Q = reshape(Q, n*size(A));