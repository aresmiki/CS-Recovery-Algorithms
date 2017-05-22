%% CS_RSL0  Robust Smoothed l0-Pseudonorm Algorithm
% 输入：y---测量信号 M X 1
%          A---恢复矩阵 M X N
%          deltaT---最小的delta 限制delta不能太小,delta越小越接近于l0范数
%          r---delta的收缩步长
%          mu---很小的更新参数，正数 2-4取值 见文章
%          L---循环次数 默认取值3
% 输出：xs---恢复的信号 N X 1
%          valL0---每次迭代的稀疏度追踪情况
% 
%  minimize ||x||_0
%  subject to ||Ax-y||_2<e_eps
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献：H. Mohimani, M. Babie-Zadeh, and C. Jutten,
% “A fast approach for overcomplete sparse decomposition based on smoothed l0-norm,”
% IEEE Trans. Signal Process., vol. 57, no. 1, pp. 289-301, Jan. 2009.
%---------------------------------------------------------------------------------------------------------------------%
%%
function [xs,valL0]=CS_RSL0(y,A,deltaT,r,mu,L,e_eps);
% r为delta缩小步长，<1。
y=y(:);
N=max(size(A));
M=min(size(A));
%% 初始化一些参数
pinv_A=pinv(A);
% pinv_A=A'*inv(A*A');
xs=pinv_A*y;  %得到一个最小二乘解作为初始值
delta=2*max(abs(xs));  %文章提出2-4倍，如果delta>4*max(abs(ys)),exp(-s.^2/(2*delta.^2))=1
k=0; %记录循环次数
valL0(1)=0;
%  maximizing F_delta(x)=sum(f_delta(x(i)))=sum(exp(-x(i).^2)/(2*delta.^2))
%  subject Ax=y;
%  用Lagrangian 推导：L(x,lambda)= F_delta(x)-lambda'*(Ax-y)
%  对x和lambda分别求导数得到KKT条件
%  x偏导数：[x(1)exp(-x(1).^2)/(2*delta.^2)),...,x(i)exp(-x(i).^2)/(2*delta.^2))-A'*lambda                (1)
%  lambda偏导数：Ax-y=0
% Ax=y的最小二乘解为：x=pinv(A)*y,相当于对优化方程
% maximizing 0.5*x*x'
% subject Ax=y
% 用Lagrangian 推导 L(x,lambda)= 0.5*x*x'-lambda'*(Ax-y)
%  x偏导数：[exp(-x(1).^2)/(2*delta.^2)),...,exp(-x(i).^2)/(2*delta.^2))-A'*lambda                (2)
%  lambda偏导数：Ax-y=0
%  (1)和(2)对比，发现delta趋近无穷大，exp((-x.^2)/(2*delta.^2))=1，两个优化问题的解几乎一样
%  delta>>max(xs),xs是方程的最小二乘解
while delta>deltaT
    k=k+1;
for i=1:L   %L次最速上升算法
    t_delta=xs.*exp(-abs(xs).^2/(2*delta^2));  %更新梯度值
    xs=xs-mu*t_delta;  %固定步长的梯度上升，获取上升后的函数值
    newy=A*xs-y;
    norm(newy,2)
    if norm(newy,2)>e_eps;   %投影
        xs=xs-pinv_A*newy; %投影到可行集上
        break;
    end

end
   valL0(k+1)=N-sum(exp(-abs(xs).^2/(2*deltaT^2)));  %要用此公式得到好的值deltaT不能太小
   delta = delta * r;  %收缩步长,逐渐缩小delta，越小越好，只要不小于deltaT
  if (abs(valL0(k+1)-valL0(k))<1e-4);
     valL0=valL0(2:end);
     break;
 end
end
valL0=valL0(2:end);




