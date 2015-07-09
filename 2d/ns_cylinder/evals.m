tic
load linear.mat;

sym = 'o';

n = 200;       % no of eigenvalues to compute
opts.p = 4*n; % size of krylov space
sigma = 20.0;

[V,D,flag] = eigs(A,M,n,sigma,opts);
assert(flag==0)
D = diag(D);
plot(real(D),imag(D),sym,'LineWidth',1.5)
grid on
D(1:10)
toc
