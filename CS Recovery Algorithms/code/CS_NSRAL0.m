%% CS_NSRAL0  Null-Space Reweigthted Approximate l0-Pseudonorm Algorithm
% 输入：y---测量信号 M X 1
%          A---恢复矩阵 M X N
%          deltaT---最小的delta 限制delta不能太小,delta越小越接近于l0范数
%          r---delta的收缩步长
%          t---很小的一个正参数，确保开始的delta在最大值上稍微大一点
%          eps---权重更新的一个较小的量eps，确保分母不为0
% 输出：yk---恢复的信号 N X 1
%          valL0---每次迭代的稀疏度追踪情况
% 
%  minimize ||x||_0
%  subject to Ax=y
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献：J. K. Pant, W.-S. Lu, and A. Antoniou, 
% “Reconstruction of sparse signals by minimizing a re-weighted approximate l0-norm in the
% null space of the measurement matrix,” IEEE Inter. Midwest Symp. on Circuits-Syst, pp. 430–433, 2010.
%---------------------------------------------------------------------------------------------------------------------%
%%
function [yk,valL0]=CS_NSRAL0(y,A,deltaT,r,t,eps);
% r为delta缩小步长，<1。
y=y(:);
N=max(size(A));
M=min(size(A));
% A维度为M X N；其零空间 V维度为 N X (N-M)   零空间的一个向量epsig为 (N-M) X 1
%% 初始化Null Space 零空间中一个向量，为通解做准备
epsig=zeros((N-M),1);   %零空间向量  
ys=A'*inv(A*A')*y; %一个最小二乘解是特解
w=ones(N,1);     %权重矩阵,开始每个权重都为1
delta=max(y)+t;  %a reasonable initial value of theta 确保开始函数是凸的
k=0;
% V=null(A);   %matlab 可以直接命令求解
[Q,R]=qr(A');   %求QR分解
V=Q(:,M+1:N);  % Q的最后N-M是矩阵A的零空间矩阵
valL0(1)=0;
 %% 构造优化问题
while (delta>deltaT)
    k=k+1;  %记录迭代次数
%  epsig=fminunc(@(epsigx) w'*(1-exp(-(ys+V*epsigx).^2./(2*delta.*delta)')),epsig);   %求最优化问题
 [epsig,val,iters]=bfgs(@(epsigx) w'*(1-exp(-(ys+V*epsigx).^2./(2*delta.*delta)')),...
                   @(epsigx) (V'*(w.*((ys+V*epsigx).*exp(-(ys+V*epsigx).^2./(2*delta.*delta)')))),epsig);
 yk=ys+V*epsig;
 w=1./(abs(yk)+eps);  %权重更新
 delta=r*delta;
 [valL0(k+1)]= ones(1,N)*(1-exp(-(yk).^2./(2*deltaT.*deltaT)'));
 if (abs(valL0(k+1)-valL0(k))<1e-4);
     valL0=valL0(2:end);
     break;
 end
end
valL0=valL0(2:end);
end




