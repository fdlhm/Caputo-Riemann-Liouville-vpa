% multND1.m
% If you use it, please cite the corresponding paper:
% Francisco de la Hoz, Peru Muniain,
% Numerical approximation of Caputo-type advection-diffusion equations via
% shifted Chebyshev polynomials, (2026).
function B=multND1(A,X)
sizeX=size(X); % size of $\mathbf X$
sizeB=[size(A,1),sizeX(2:end)]; % size of $\mathbf B$
B=reshape(A*reshape(X,sizeX(1),[]),sizeB); % $\mathbf B$