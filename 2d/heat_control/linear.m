% State-space form
% E dx/dt = A x + Bu
clear all
load linear.mat

[nu,nl] = size(N);

A = [A, N;
      N', zeros(nl,nl)];
E = [M, zeros(nu,nl);
     zeros(nl,nu+nl)];
B = [zeros(nu,nl);
     -Mb];

n = 10;
opts.p = 4*n;

fprintf(1,'Eigenvalues of A,E\n')
[V,D,flag] = eigs(A, E, n, 'SM', opts);
D = diag(D)/pi^2;
fprintf(1,'Eigenvalues/(pi*pi) =\n')
fprintf(1,'%f\n',D)

% Adjoint
fprintf(1,'Eigenvalues of transpose\n')
[V,D,flag] = eigs(A', E', n, 'SM', opts);
D = diag(D)/pi^2;
fprintf(1,'Eigenvalues/(pi*pi) =\n')
fprintf(1,'%f\n',D)

% Check Hautus if first eigenvalue is unstable
h = B' * V(:,1);
fprintf(1,'Norm of hautus = %e\n',norm(h))
