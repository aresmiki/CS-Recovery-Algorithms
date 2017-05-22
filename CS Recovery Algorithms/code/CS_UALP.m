%% CS_UALP  Minimization of Approximate Lp Pseudonorm Using a Quasi-Newton Algorithm
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
%  subject to Ax=y
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献：Pant J K, Lu W S, Antoniou A. 
% Unconstrained regularized Lp -norm based algorithm for the reconstruction of sparse signals[C]
% IEEE International Symposium on Circuits and Systems. IEEE, 2011:1740-1743.
%---------------------------------------------------------------------------------------------------------------------%
%% UALP是Pant J K在博士论文中的定义，URLP是IEEE论文中的定义，算法一样
function yk=CS_UALP(y,A,p)
y=y(:);
N=max(size(A));
M=min(size(A));
% A维度为M X N；其零空间 V维度为 N X (N-M)   零空间的一个向量epsig为 (N-M) X 1
%% 初始化Null Space 零空间中一个向量，为通解做准备
epsig=zeros((N-M),1);   %零空间向量  
ys=A'*inv(A*A')*y; %一个最小二乘解是特解
% k=0;
% V=null(A);   %matlab 可以直接命令求解
[Q,R]=qr(A');   %求QR分解
V=Q(:,M+1:N);  % Q的最后N-M是矩阵A的零空间矩阵
eps1=sqrt(1-p)*max(abs(ys));
epsT=1e-5;
T=9;
belta=log(eps1/epsT)/(T-1);
for i=2:(T-1)
    eps(i)=exp(-belta*i);
end
eps=[eps1,eps,epsT];
w=ones(N,1);     %权重矩阵,每个权重都为1,为了求和
for k=1:T
      epsig=fminunc(@(epsigx) w'*((ys+V*epsigx).^2+eps(k).^2).^(p/2),epsig);   %求最优化问题
%     [epsig,val,iters]=bfgs(@(epsigx) w'*((ys+V*epsigx).^2+eps(k).^2).^(p/2),...
%                    @(epsigx) p*V'*((((ys+V*epsigx).^2+eps(k).^2).^(p/2-1)).*(ys+V*epsigx)),epsig);             
end
yk=ys+V*epsig;



%%  BFGS中的线性搜索
%     function [alpha1]=LSBFP(epsigk,xs,V,p,dk,deltaT,epsk)
%         lalpha1=0;deltaA=deltaT+1;
%         while(deltaA>deltaT)
%             gy=(xs+V*epsigk+alpha1*V*dk).^2+epsk*epsk;
%             Geps=-sum(((xs+V*epsigk).*(V*dk).*gy.^(p/2-1)))/sum(V*dk.*gy.^(p/2-1));
%             deltaA=Geps-alpha1;
%             alpha1=Geps;
%         end
%         
%     end

end

