%% Iterative Shrinkage-Thresholding Algorithm    (ISTA)
%----------------------- A Proximal-Gradient Algorithm Method-------------------------%
%  CS_ISTA  Algorithm (迭代收缩阈值算法 ISTA)   
%  输入：y---测量信号  M X 1
%           A---恢复矩阵  M X N
%           lambda---正则化参数     
%           iter---最大迭代次数
% 输出 ：xhk---估计的稀疏向量 N X 1
%            err---每次迭代的误差
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月26日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献：I. Daubechies, M. Defrise, and C. D. Mol,
% “An iterative thresholding algorithm for linear inverse problems with a sparsity constraint,” 
% Comm. Pure Appl. Math., vol. 57, pp. 1413-1457, 2004.
%------------------------------------------------------------------------------------------%
%  minimize ||A*x-y||^2 + lambda*||x||_1     (1)
%  ||A*x-y||^2是二范数平方
%  f(x)= ||A*x-y||^2=x'A'Ax-2x'A'y+y'y
%  一阶导数为：▽f(x)= 2A'Ax-2A'y=2A'(Ax-y)

% (1)的最小化问题转换为 x_k=argmin{(1/(2*t_k))||x-(x_k-t_k▽f(x_k-1))||^2+lambda||x||_1}
% 等价为x_k=argmin{(1/(2*t_k))||x-c_k||^2+lambda||x||_1}
% 其中c_k=x_k-t_k▽f(x_k-1)=x_k-1-2t_kA'(Ax_k-1-y)

% t_k为求解步长
%%
%t_k=1/L(f)  Lipschitz常数的倒数
% 算法实现,c_k,t_k,lambda,alp=lambda*t_k; t_k=[0,1/||A'A||]
% x_k=max(abs(c_k)-alp,0).*sign(c_k)
%%

% F(x_k)-F(x*)<=L(f)||x_0-x*||^2/2*k
%收敛速度1/k

%%
function [xhk,err]=CS_ISTA(y,A,lambda,iter)
if nargin<4;
    iter=1e6;  %足够大，直到能达到精度1e-6
end
y=y(:);
N=max(size(A));
M=min(size(A));
% tmax=1/max(eig(A*A'));   %求解Lipschitz常数   1/||A'A||
%-------------------------------------------%
% 用power方法求最大特征值
%
Mat=A*A';
x=randn(M,1);
for i=1:5;  %一般3-5次就可以了
    x=x/norm(x,2); %归一化
    x=Mat*x;
end
x=x/norm(x,2);
Lf=x'*Mat*x;
tmax=1/Lf;
%-------------------------------------------%

tk=0.95*tmax;  %上界稍微变小一点
xhk=zeros(N,1);  %初始化存储向量，开始迭代的值从0开始 相当于x_k；
% xhk=A'*y;       %可以估计一个最小二乘值，从该值开始迭代
alp=lambda*tk;  %计算alp=lambda*t_k
err(1)=0;
for i=1:iter;
     ck=xhk-2*tk*A'*(A*xhk-y);  %计算c_k
     xhk=(max(abs(ck)-alp,0)).*sign(ck);  %更新x_k,进入下次迭代
     err(i+1)=norm(A*xhk-y,2);
     if abs(err(i+1)-err(i))<1e-6;  %跳出循环
         err=err(2:end);
         break;
     end
end
err=err(2:end);       
end