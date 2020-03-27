%问题19的解析接
%首先对矩阵做特征值分解
%将分解出来的对角阵设置阈值即可求解出结果。
%由于eig，命令是将特征值按照升序排列的，无需对Dga改变顺序
function Z=step2_closedform(N,x)
gamma=0.1;
y=x';
y(1)=x(1);
T=toeplitz(x,y);
[U,D]=eig(T);
Dga=zeros(N,N);
for i=1:N
    if(gamma-D(i,i)>0)
        Dga(i,i)=gamma-D(i,i);
    else
        Dga(i,i)=0;
    end
end
Z=U*Dga*U';
%求出的Z即为问题19的全局最优解 这便是完整过程
end