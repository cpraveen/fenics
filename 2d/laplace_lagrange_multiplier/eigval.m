load linear.mat

[nu,nl] = size(N);

A1 = [A, N;
      N', zeros(nl,nl)];
M1 = [M, zeros(nu,nl);
      zeros(nl,nu+nl)];

n = 10;
opts.p = 4*n;
[V,D,flag] = eigs(A1, M1, n, 'SM', opts);
D = diag(D)/pi^2;
fprintf(1,'Eigenvalues/(pi*pi) =\n')
fprintf(1,'%f\n',D)
