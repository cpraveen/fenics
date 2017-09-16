function r = is_diag(A)

tmp = A - diag(diag(A),0);
tmp = max(max(abs(tmp)));
if tmp<1.0e-10
   r = 1;
else
   r = 0;
end
