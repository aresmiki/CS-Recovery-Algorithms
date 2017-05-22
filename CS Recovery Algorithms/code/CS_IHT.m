%% CS_IHT  Iterative Hard Thresholding algorithms for compressive sensing
% 输入：y---测量信号 M X 1
%          A---恢复矩阵 M X N
%          K---信号的稀疏度
% 输出：x0---恢复的信号 N X 1
%    
% 
%  minimize ||x||_1
%  subject to ||Ax-y||_2<eps
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献：Blumensath T, Davies M E. 
% Iterative hard thresholding for compressed sensing[J]. 
% Applied & Computational Harmonic Analysis, 2009, 27(3):265-274.
% Blumensath T, Davies M E.
% Iterative Thresholding for Sparse Approximations[J].
% Journal of Fourier Analysis and Applications, 2008, 14(5):629-654.
%---------------------------------------------------------------------------------------------------------------------%

function x0=CS_IHT(y,A,K)

M=min(size(A));
N=max(size(A));
if nargin<3;
    K=floor(M/4);        %最少迭代次数，一般等于稀疏度数，经验公式除以4，为保证重构精度，iter可以选大一点
end
x0=zeros(N,1);         % 初始化解空间向量
u=0.5;                       % 影响系数

for times=1:M
    
    x1=A'*(y-A*x0);
    
    x2=x0+u*x1;
    
    [val,pos]=sort(abs(x2),'descend');  % 降序排列
    
    x2(pos(K+1:end))=0;   % 保留最大的前iters个数的数据，iters为稀疏度

    x0=x2;          % 更新值到下一步循环

end
