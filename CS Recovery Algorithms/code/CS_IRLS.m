%% CS_IRLS  Iteratively reweighted algorithms for compressive sensing
% 输入：y---测量信号 M X 1
%          A---恢复矩阵 M X N
%          p---伪范数 1--0之间
%          deltaT---最小的delta 限制delta不能太小,delta越小越接近于l0范数
%          r---delta的收缩步长
%          mu---很小的更新参数，正数 2-4取值 见文章
%          L---循环次数 默认取值3
% 输出：ys---恢复的信号 N X 1
%    
% 
%  minimize ||x||_p
%  subject to Ax=y
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月30日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献： Chartrand and W. Yin,
% “Iteratively Reweighted Algorithms for Compressed Sensing,” 2008.
%---------------------------------------------------------------------------------------------------------------------%
function ys=CS_IRLS(y,A,p,epsT,r)
N=max(size(A));
M=min(size(A));
y=y(:);

eps=1;
ys=inv(A'*A)*A'*y;  %估计一个最小二乘解
% ys=pinv(A)*y;

while (eps>epsT)
    
    w=(ys.^2+eps).^(p/2-1);  %更新权重
%     w=(abs(ys)+eps).^(p-1);  %更新权重
    Q=diag(1./w);
    yx=Q*A'*inv(A*Q*A')*y;  %更新
    if(norm(yx-ys,2) < sqrt(eps)*r.^2)
        eps=eps*r;  %更新eps  r倍更新
    end
    ys=yx;
    
end

end
