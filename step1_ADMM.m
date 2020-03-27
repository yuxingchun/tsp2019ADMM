function [u,x,P1,P2,B1,B2]=step1_ADMM(y,Z,u,x,P1,P2,B1,B2,eta,N,rho)
Tu=toeplitz(u);
irho=1/rho;

%updata u
ResTu=-1.*irho.*mapM2v(Z,N)+irho.*mapM2v(B2(1:N,1:N),N)+mapM2v(P2(1:N,1:N),N);
u=imapM2v(ResTu,N);

%updata x
x=0.5.*y+0.5.*irho.*(B1(1:N,N+1)+B2(1:N,N+1));+0.5.*(P1(1:N,N+1)+P2(1:N,N+1));

%updata P1
P1=[eta.*eye(N),x-y;(x-y)',eta]-irho.*B1;

%updata P2
P2=[Tu,x;x',trace(Tu)/N]-irho.*B2;

%updata Q1 Q2 B1 B2
Q1=[eta.*eye(N),x-y;(x-y)',eta];
Q2=[Tu,x;x',trace(Tu)/N];
B1=B1+rho.*(P1-Q1);
B2=B2+rho.*(P2-Q2);

end

function res=mapM2v(W,M)
v = zeros(M,1);
v(1) = trace(W);
for j = 1:M-1
    v(j+1) = 2 * sum(diag(W,j));
end
res = v;
end

function ires=imapM2v(v,M)

weight = 1 ./ [M; 2*(M-1:-1:1)'];
ires=weight.*v;
end