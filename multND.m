% multND.m
% If you use it, please cite the corresponding paper:
% Carlota M. Cuesta, Francisco de la Hoz,
% A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
% arXiv:2412.15840, (2024).
%
% Multiply a matrix and an N-dimensional array along the jth dimension
function B=multND(A,X,j) % compute $\mathbf B=\mathbf A\square_j\mathbf X$
if isscalar(A) % test whether A is a scalar
    B=A*X; % trivial case
else
    sizeX=size(X); % size of X
    permutation=1:length(sizeX); % length(sizeX) is the number of dimensions N
    permutation([1 j])=permutation([j 1]); % permutation of the dimensions
    X_permuted=permute(X,permutation); % swap the first and jth dimensions
    X_reshaped=reshape(X_permuted,sizeX(j),[]); % reshape X
    B_reshaped=A*X_reshaped; % compute the product
    B_permuted=reshape(B_reshaped,size(X_permuted)); % reshape the product
    B=permute(B_permuted,permutation); % swap the first and jth dimensions
end
