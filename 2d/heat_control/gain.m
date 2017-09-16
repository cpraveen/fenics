% State-space form
% E dx/dt = A x + Bu
clear all
load linear.mat

mode = 'lqr'

% nt = no. of temp dofs
% nl = no. of lagrange mult
%    = no. of boundary control variables
[nt,nl] = size(A12);

A = [A11, A12;
     A12', sparse(nl,nl)];
E = [E11, sparse(nt,nl);
     sparse(nl,nt+nl)];
B = [sparse(nt,nl);
     -Mb];

ne = 10; % No. of eigenvalues to compute
opts.p = 4*ne;

fprintf(1,'Eigenvalues of A,E\n')
[V,D,flag] = eigs(A, E, ne, 'SM', opts);
D = diag(D);
fprintf(1,'Eigenvalues/(pi*pi) =\n')
fprintf(1,'%f\n',D/(pi^2))

% find unstable eig
iu = find(real(D) > 0);
nu = length(iu);
fprintf(1, 'Number of unstable eigenvalues of A = %d\n', nu)

% CHECK: Check Hautus if first eigenvalue is unstable
h = B' * V(:,1);
fprintf(1,'Norm of hautus = %e\n',norm(h))

% Get unstable eigenvectors
Vt = V(:,1);
Zt = V(:,1);

% CHECK: p must be diagonal
disp('Following must be a diagonal matrix. Is it ?')
p = Vt.' * E * Zt
% Check p is diagonal and diagonal entries are non-zero
assert(min(abs(diag(p))) > 0.0);
assert(is_diag(p)==1)
p = diag(p);

% Eigenvector corresponding to lagrange multiplier
Zp = Zt(nt+1:end,:);
% Eigenvectors for temperature
Vy = Vt(1:nt,:);
Zy = Zt(1:nt,:);
B1 = sparse(nt,nl);
B2 = Mb;

% CHECK: orthonormality
disp('Is this identity matrix ?')
p = Vy.' * E11 * Zy
% Check diagonal entries are 1
assert(max(abs(diag(p)-1.0)) < 1.0e-10);
assert(is_diag(p)==1)

N  = [E11, A12; A12' sparse(nl,nl)];
RHS= [sparse(nt,nl); B2];
Z1 = N\RHS;
B12= B1 + A11*Z1(1:nt,:);

% Project to unstable subspace
Au = Zy' * A11 * Vy;
Bu = Zy' * B12;

if mode == 'min'
   % minimal norm control
   Ru = eye(nl);
   Qu = zeros(size(Au));
   fprintf(1,'Minimal norm feedback\n')
else
   % LQR problem
   Ru = eye(nl);
   Qu = Vy' * E11 * Vy;
   fprintf(1,'LQR feedback\n')
end

[Pu,L,G]=care(Au,Bu,Qu,Ru);
disp('Eigenvalues of projected system with feedback')
L
disp('Eigenvalues of Pu')
ePu = eig(Pu)
max_ePu = max(ePu);
% save to file
fid=fopen('maxeig.dat','w');
fprintf(fid,'%24.14e',max_ePu);
fclose(fid);

B = sparse([B1; -B2]);
E11 = sparse(E11);
Z = sparse([Zy; Zp]);
Zy = sparse(Zy);
Pu = sparse(Pu);
Kt = Ru \ ((B' * Z) * Pu * (Zy' * E11));

% CHECK: are we stable now
S  = [Kt, sparse(nl,nl)];
A  = sparse([A11, A12; A12', sparse(nl,nl)]);
M  = sparse([E11, sparse(nt,nl); sparse(nl,nt+nl)]);
[V,D,flag] = eigs(A-B*S,M,ne,'SM',opts);
assert(flag==0)
disp('Eigenvalues of full system with feedback')
diag(D)/pi^2

K = full(Kt);
save('gain.mat','K');
