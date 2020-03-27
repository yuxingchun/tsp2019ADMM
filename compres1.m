function [res_prim,res_dual,tol_prim,tol_dual,U]=compres1(Z,u,x,P2,B2,eta,N,rho,U)
tol_abs = 1e-4;
tol_rel = 1e-5;
U_last=U;
psdmat = [toeplitz(u), x; x', (N^-1)*trace(toeplitz(u))];
U_temp = psdmat + 1/rho .* B2;
[p, q] = eig((U_temp + U_temp') / 2); 
%找到实部大于0的指标
dq = real(diag(q));
idxpos = (dq > 0);
U = p(:,idxpos) * diag(dq(idxpos)) * p(:,idxpos)';
res_prim = psdmat - U;
U_dif = U - U_last;
res_dual = rho * [U_dif(N+1,N+1); 2*U_dif(1:N,N+1); mapM2v(U_dif(1:N,1:N),N)];

tol_prim = sqrt(2)*(N+1)*tol_abs + tol_rel*max(norm(psdmat,'fro'), norm(U,'fro'));
tol_dual = sqrt(2*(N+1))*tol_abs + tol_rel*norm([B2(N+1,N+1); 2*B2(1:N,N+1); mapM2v(B2(1:N,1:N),N)]);
end

function res=mapM2v(W,M)
v = zeros(M,1);
v(1) = trace(W);
for j = 1:M-1
    v(j+1) = 2 * sum(diag(W,j));
end
res = v;
end