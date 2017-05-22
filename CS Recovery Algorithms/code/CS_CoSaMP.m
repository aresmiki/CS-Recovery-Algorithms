%% CS_CoSaMP  Algorithm 
%-----------------------------------------------------------------------------------------%
%  CS_CoSaMP  Algorithm (压缩采样正交匹配追踪法 Orthogonal Matching Pursuit)   
%  输入：y---测量信号  M X 1
%           A---恢复矩阵  M X N
%           K---迭代次数
% 输出 ：theta---估计的稀疏向量 N X 1
%            erro_res---每次迭代的误差
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年04月26日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%  参考文献1：Needell D，Tropp J A 
%  CoSaMP：Iterative signal recovery from incomplete and inaccurate samples[J]．
% Applied and Computation Harmonic Analysis，2009，26：301-321.
% 参考文献2：D.Needell, J.A. Tropp．
% CoSaMP: Iterative signal recoveryfrom incomplete and inaccurate samples[J]. 
% Communications of theACM，2010，53(12)：93-100.
%------------------------------------------------------------------------------------------%

%%
function [ theta,erro_res ] = CS_CoSaMP( y,A,K )
    [m,n] = size(y);
    if m<n
        y = y'; 
    end
    [M,N] = size(A); %传感矩阵A为M*N矩阵
    theta = zeros(N,1); %用来存储恢复的theta(列向量)
    pos_num = []; %用来迭代过程中存储A被选择的列序号
    res = y; %初始化残差(residual)为y
    for kk=1:K %最多迭代K次
        %(1) Identification
        product = A'*res; %传感矩阵A各列与残差的内积
        [val,pos]=sort(abs(product),'descend');
        Js = pos(1:2*K); %选出内积值最大的2K列
        %(2) Support Merger
        Is = union(pos_num,Js); %Pos_theta与Js并集
        %(3) Estimation
        %At的行数要大于列数，此为最小二乘的基础(列线性无关)
        if length(Is)<=M
            At = A(:,Is); %将A的这几列组成矩阵At
        else %At的列数大于行数，列必为线性相关的,At'*At将不可逆
            if kk == 1
                theta_ls = 0;
            end
            break; %跳出for循环
        end
        %y=At*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y; %最小二乘解
        %(4) Pruning
        [val,pos]=sort(abs(theta_ls),'descend');
        %(5) Sample Update
        pos_num = Is(pos(1:K));
        theta_ls = theta_ls(pos(1:K));
        %At(:,pos(1:K))*theta_ls是y在At(:,pos(1:K))列空间上的正交投影
        res = y - At(:,pos(1:K))*theta_ls; %更新残差 
        erro_res(kk)=norm(res,2);
        if norm(res)<1e-6 %Repeat the steps until r=0
            break; %跳出for循环
        end
    end
    theta(pos_num)=theta_ls; %恢复出的theta
end