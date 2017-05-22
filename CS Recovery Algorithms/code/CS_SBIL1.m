%%  CS_SBIL1  L1_SplitBregmanIteration
% 输入：f---测量信号 M X 1
%          A---恢复矩阵 M X N
%          mu---步长
%          lambda---规则化因子
%          Niter---最大迭代次数
% 输出：u---恢复的信号 N X 1
%    
% 
%  minimize ||x||_1
%  subject to Ax-y=0
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
% 参考文献：Yin W, Osher S, Goldfarb D, et al. 
% Bregman Iterative Algorithms for L1 Minimization with Applications to Compressed Sensing[J].
% Siam Journal on Imaging Sciences, 2008, 1(1):143-168.
%---------------------------------------------------------------------------------------------------------------------%


function u=CS_SBIL1(f,A,mu,lambda,Niter)

N=max(size(A));

d=zeros(N,1);
b=zeros(N,1);
u=zeros(N,1);

Z=zeros(N,1);
Ft=mu*A'*f;
IV=inv(mu*(A'*A)+lambda*eye(N));

err=norm(f,2);
tol=1e-3*err;
K=0;
while ((err>tol) && (K<Niter)),
    K=K+1;
    up=u;
    u=IV*(Ft+lambda*(d-b));
    tmp=u+b;
    d=sign(tmp).*max(Z,abs(tmp)-1/lambda);
    b=tmp-d;
    err=norm(u-up,2);
end