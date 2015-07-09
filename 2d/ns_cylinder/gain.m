clear all
load linear.mat
load freeinds.txt
load pinds.txt
who

% Set some parameters
mode = 'min'; % 'min' or 'lqr'
ne = 50;      % how many eigenvalues to compute
nes = 2;      % number of eigenvalues to shift
% With shift=1, there are 6 complex unstable eigenvalues
shift = 0.0;

nc = size(B,2);   % number of control variables

% number of Lanczos vectors
opts.p = 4*ne;
sigma  = 10.0;

% Compute eigenvalues,vectors of (A,M)
[Vt,D1,flag] = eigs(A,M,ne,sigma,opts);
assert(flag==0)
disp('Eigenvalues of A')
D1=diag(D1)

% find unstable eig
iu = find(real(D1) > 0);
nu = length(iu);
fprintf(1, 'Number of unstable eigenvalues of A = %d\n', nu)

% select nes eigenvalues with largest real part
% NOTE: Make sure they appear in conjugate pairs
[d,ii]=sort(real(D1), 'descend');
D1 = D1(ii(1:nes))
% with shift all the D1 must be unstable
assert(min(real(D1+shift)) > 0.0)
Vt = Vt(:,ii(1:nes));

% Compute eigenvalues,vectors of (A^T,M^T)
[Zt,D2,flag] = eigs(A',M',ne,sigma,opts);
assert(flag==0)
disp('Eigenvalues of A^T')
D2=diag(D2)

% find unstable eig
iu = find(real(D2) > 0);
nu2= length(iu);
fprintf(1, 'Number of unstable eigenvalues of A^T= %d\n', nu2)
assert(nu == nu2)

% select nes eigenvalues with largest real part
[d,ii]=sort(real(D2), 'descend');
D2 = D2(ii(1:nes))
% with shift all the D2 must be unstable
assert(min(real(D2+shift)) > 0.0)
Zt = Zt(:,ii(1:nes));

% NOTE: check that eigenvalues are in same order
for j=1:nes
   assert(abs(D1(j) - D2(j)) < 1.0e-10)
   % check all eigenvalues are complex
   assert(abs(imag(D1(j))) > eps)
   Vt(:,j) = Vt(:,j) / max(abs(Vt(:,j)));
end

% make Vt and Zt orthonormal
% p must be diagonal
disp('Following must be a diagonal matrix. Is it ?')
p = Vt.' * M * Zt
% Check p is diagonal and diagonal entries are non-zero
assert(min(abs(diag(p))) > 0.0);
assert(is_diag(p)==1)
p = diag(p);

% normalize
for j=1:nes
   Zt(:,j) = Zt(:,j) / p(j);
end

% freeinds, pinds are indices inside fenics
% We have to shift by one since python indexing starts at 0 but matlab 
% starts at 1
freeinds = freeinds + 1;
pinds    = pinds + 1;
% get matlab indices of velocity+temperature
[tmp,vinds] = setdiff(freeinds, pinds, 'stable');
% get matlab indices of pressure
pinds = setdiff(1:length(freeinds), vinds, 'stable');

% eigenvector component for velocity+temperature
Vty = Vt(vinds,:);  % eigenvectors of (A,M)
Zty = Zt(vinds,:);  % eigenvectors of (A',M')
Ztp = Zt(pinds,:);

E11 = M(vinds,vinds);
A11 = A(vinds,vinds);
A12 = A(vinds,pinds);
B1  = B(vinds,:);
B2  =-B(pinds,:);

% check orthonormality
disp('Is this identity matrix ?')
p = Vty.' * E11 * Zty
% Check diagonal entries are 1
assert(max(abs(diag(p)-1.0)) < 1.0e-10);
assert(is_diag(p)==1)

U = (1/sqrt(2)) * [1,   1; ...
                   1i, -1i];
%U = blkdiag(U,U,U);
assert(size(U,1) == nes)

Vy = Vty * U';
Zy = Zty * U.';
Zp = Ztp * U.';

disp('Vy and Zy must be real')
assert(max(max(abs(imag(Vy)))) < 1.0e-13)
assert(max(max(abs(imag(Vy)))) < 1.0e-13)

% Vy and Zy must be real, making sure imaginary part is close to zero
Vy = real(Vy);
Zy = real(Zy);
Zp = real(Zp);

disp('Is this identity matrix ?')
p = Vy.' * E11 * Zy
% Check diagonal entries are 1
assert(max(abs(diag(p)-1.0)) < 1.0e-10);
assert(is_diag(p)==1)

% Compute B12
np = length(pinds);
ny = length(vinds);
N  = [E11, A12; A12' sparse(np,np)];
RHS= [sparse(ny,nc); B2];
Z1 = N\RHS;
B12= B1 + A11*Z1(1:ny,:);

% Project to unstable subspace
Au = Zy' * A11 * Vy;
Au = Au + shift*eye(size(Au));
Bu = Zy' * B12;

if mode == 'min'
   % minimal norm control
   Ru = eye(nc);
   Qu = zeros(size(Au));
   fprintf(1,'Minimal norm feedback\n')
else
   % LQR problem
   Ru = diag([0.05, 0.01, 0.05]);
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
% uncomment following if you are running best_loc.py
%exit

B = sparse([B1; -B2]);
E11 = sparse(E11);
Z = sparse([Zy; Zp]);
Zy = sparse(Zy);
Pu = sparse(Pu);
Kt = Ru \ ((B' * Z) * Pu * (Zy' * E11));
S  = [Kt, sparse(nc,np)];
A  = sparse([A11, A12; A12', sparse(np,np)]);
M  = sparse([E11, sparse(ny,np); sparse(np,ny+np)]);
[V,D,flag] = eigs(A-B*S,M,ne,sigma,opts);
assert(flag==0)
disp('Eigenvalues of full system with feedback')
diag(D)

Kt = full(Kt);
save('gain.mat','Kt')
save('state.mat','M','A','B','S','Z')
