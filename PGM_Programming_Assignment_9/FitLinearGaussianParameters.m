function [Beta sigma] = FitLinearGaussianParameters(X, U)

% Estimate parameters of the linear Gaussian model:
% X|U ~ N(Beta(1)*U(1) + ... + Beta(n)*U(n) + Beta(n+1), sigma^2);

% Note that Matlab/Octave index from 1, we can't write Beta(0).
% So Beta(n+1) is essentially Beta(0) in the text book.

% X: (M x 1), the child variable, M examples
% U: (M x N), N parent variables, M examples
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

M = size(U,1);
N = size(U,2);

Beta = zeros(N+1,1);
sigma = 1;

% collect expectations and solve the linear system
% A = [ E[U(1)],      E[U(2)],      ... , E[U(n)],      1     ; 
%       E[U(1)*U(1)], E[U(2)*U(1)], ... , E[U(n)*U(1)], E[U(1)];
%       ...         , ...         , ... , ...         , ...   ;
%       E[U(1)*U(n)], E[U(2)*U(n)], ... , E[U(n)*U(n)], E[U(n)] ]

% construct A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N+1,N+1);
for c=1:N
    A(1,c) = mean(U(:,c));
end
A(1, N+1) = 1;
for r=1:N
    A(r+1,N+1) = mean(U(:,r));
end

for r=1:N
    for c=1:N
        A(r+1,c) = mean(U(:,c) .* U(:,r));
    end
end

% B = [ E[X]; E[X*U(1)]; ... ; E[X*U(n)] ]

% construct B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(1, N+1);
B(1,1) = mean(X);
for c=1:N
    B(c+1, 1) = mean(X.*U(:,c));
end
% solve A*Beta = B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta = A \ B;
Beta = Beta(:,1);
% then compute sigma according to eq. (11) in PA description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cov_xx = cov(X, 1);
%cov_xx = mean(X .* X) - mean(X) * mean(X);
t = 0.0;
for i=1:N
    for j=1:N
        %cov_uij = cov(U(:,i), U(:,j), 1);
        cov_uij = mean(U(:,i) .* U(:,j)) - mean(U(:,i)) * mean(U(:,j));
        temp = Beta(i)*Beta(j)*cov_uij;
        t = t + temp;
    end
end

sigma = sqrt(cov_xx - t);

