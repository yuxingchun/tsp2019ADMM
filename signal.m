
function y0=signal(N)
%信号数量K
K=4;
%频率 加噪声
f = [0.103; 0.115; 0.5;0.8] + unifrnd(-.001,.001,K,1);
%幅度
am = [2; 2; 1;0.5];
%幅度加入相位
s = am .* exp(1i*2*pi*rand(K,1)); 
Omega=0:N-1;
Omega=Omega';
A = exp(1i*2*pi*kron(Omega,f')); 
y0 = A * s; 
end
