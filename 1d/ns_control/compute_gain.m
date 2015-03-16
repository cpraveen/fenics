load linear.mat

M = full(M);
A = full(A);
B = [B_ul, B_ur];

e = eigs(A,M,50,'LR');

ii = find(real(e)>0);
e(ii)

figure(1)
plot(real(e), imag(e), 'o')
grid on
hold on

Q = zeros(size(A));
R = eye(2);
S = [];

[X,L,G] = care(A,B,Q,R,S,M);

plot(real(L),imag(L),'*')
axis([-50 1 -200 200])

save('gain.mat','G')
