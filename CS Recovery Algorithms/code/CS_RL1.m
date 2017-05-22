%% CS_RL1 Reweighted L1 Minimization Algorithm
% 输入：y---测量信号 M X 1
%          A---恢复矩阵 M X N
%          iter---最大迭代次数  至少大于2，小于2就相当于L1最小化（没有加权）

% 输出：ys---恢复的信号 N X 1
%    
% 
%  minimize W||x||_1
%  subject to Ax=y
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献： Candès E J, Wakin M B, Boyd S P. 
% Enhancing sparsity by reweighted L1 minimization.[J]. 
% Journal of Fourier Analysis & Applications, 2007, 14(5):877-905.
%---------------------------------------------------------------------------------------------------------------------%
function xh=CS_RL1(y,A,iter)
N=max(size(A));
M=min(size(A));
Li=(M/(4*log(N/M)));
y=y(:);
W=ones(N,1);   %初始化权重向量
QW=diag(W);    %初始化权重矩阵
% delta=0.01;
for i=1:iter
    QWt=inv(QW);
    At=A*QWt;
    x0=At'*y;  %最小二乘解估计一个初始值
    xh=l1eq_pd(x0,At,[],y,1e-3);
    delta=max(norm(xh,Li),1e-3) ;%动态更新 自适应更新delta
%     delta should be set slightly smaller than the expected nonzero magnitudes of x0. 
%     In general, the recovery process tends to be reasonably robust to the choice of delta.--原文中的话
    xh=QWt*xh;
    QW=diag(1./(abs(xh)+delta));
end

end
