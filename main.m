%输入分辨率N^-1
N=30;
maxit=2000;
SNR=5;
gamma=0.1;
%得到无噪声信号
y0=signal(N);
%加入噪声
sigma2=1/SNR;
eta=sqrt((N+2*sqrt(N))*sigma2);
%sigma = min(am).^2 / 10^(SNR/10);
y = y0 + sqrt(sigma2/2)*(randn(N,1)+1i*randn(N,1));
%初始化
Z=gamma*eye(N);
u=zeros(N,1);
x=zeros(N,1);
Tu=toeplitz(u);
P1=[eta.*eye(N),x-y;(x-y)',eta];
P2=[Tu,x;x',trace(Tu)/N];
B1=zeros(N+1);
B2=zeros(N+1);
mu=10;
tau=2;
U=zeros(N+1);
U2=zeros(N+1);
rho=2;
for i=1:maxit
%objvalue1=trace(Z*toeplitz(u));
[u,x,P1,P2,B1,B2]=step1_ADMM(y,Z,u,x,P1,P2,B1,B2,eta,N,rho);
Z=step2_closedform(N,x);
%objvalue2=trace(Z*toeplitz(u));
%终止条件 
[res_prim,res_dual,tol_prim,tol_dual,U]=compres1(Z,u,x,P2,B2,eta,N,rho,U);
err_prim = norm(res_prim, 'fro');
err_dual = norm(res_dual);
if err_prim < tol_prim && err_dual < tol_dual
        break;
end

if err_prim > mu * err_dual
        rho = tau * rho;
    elseif err_prim < err_dual / mu
        rho = rho / tau;
end
end
[freq, amp, sigma_in_u] = VanDec(u);
