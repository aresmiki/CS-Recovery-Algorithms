

%%
%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(l1-MAGIC和l1_ls解l1问题)     信号本身就是稀疏的，
%  不需要稀疏矩阵，恢复矩阵A是单位正交矩阵，用OMP方法求解l1问题.
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--西南交通大学牵引动力国家重点实验室 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月27日
%---------------------------------------------------------------------------------%
clc;clear all;close all;
%% 1. 产生稀疏的信号
N=1024;
K=50;
x=zeros(N,1);
rand('state',8)
q=randperm(N); %随机排列1到N的整数
randn('state',10)
x(q(1:K))=randn(K,1); %将K个随机数随机放到x中
t=0:N-1;
%% 2. 构造感知矩阵
M=2*ceil(K*log(N/K));
Phi=randn(M,N);  %高斯矩阵作为感知矩阵
Phi=orth(Phi')';  %正交化

%% 3. 量测信号
y=Phi*x;
A=Phi;%恢复矩阵,稀疏化矩阵为单位矩阵，因为信号本身就是稀疏的，不需要做任何稀疏变换
%% 4. 重构信号 l1最小化   Using  l1-MAGIC  该方法恢复较好 模型匹配min_x ||x||_1  s.t.  Ax = b
x0=A'*y;  %最小二乘解估计一个初始值
xh1=l1eq_pd(x0,A,[],y,1e-3);

%% 5. 重构信号l1最小化   Using l1_ls 该方法稍微差点 模型为 minimize ||A*x-y||^2 + lambda*sum|x_i|,
lambda  = 0.01; % 正则化参数
rel_tol = 1e-3; % 目标相对对偶间隙
quiet=1;   %不输出中间结果
[xh2,status]=l1_ls(A,y,lambda,rel_tol,quiet);

%% 6. 恢复信号和原始信号比较 
figure
plot(t,xh1,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-MAGIC恢复信号','原始信号')

figure
plot(t,xh2,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-ls恢复信号','原始信号')

%% 7. 用正交匹配追踪的方法计算l1最优化问题
[ xh,erro_rn ] = CS_OMP( y,A,2*K );
figure
plot(erro_rn,'-*')
legend('OMP正交匹配追踪误差')
figure
plot(t,xh,'ko',t,x,'r.')
xlim([0,t(end)])
legend('OMP恢复信号','原始信号')

%% 8. 用压缩采样匹配追踪(CoSaMP)的方法计算l1最优化问题
%Needell D, Tropp J A. CoSaMP: Iterative signal recovery from incomplete 
% and inaccurate samples [J]. Applied & Computational Harmonic Analysis, 2008, 26(3):301-321.
% 一次性选择2*K个较大的基，每次循环不断删除和补进一定数目的基，最快达到给定精度

A=Phi;    %恢复矩阵

%% CoSaMP
 [theta,erro_rnn]=CS_CoSaMP( y,A,K );
 figure
plot(erro_rnn,'-*')
legend('CoSaMP压缩采样追踪误差')
 %%    
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('CoSaMP恢复信号','原始信号')

%% CS_NSRAL0
A=Phi;    %恢复矩阵 
deltaT=1e-3;
r=1/3;
te=0.01;
eps=0.09;
[theta,Spare_L0]=CS_NSRAL0(y,A,deltaT,r,te,eps);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('近似L0范数追踪结果','真实稀疏度')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('NSRAL0恢复信号','原始信号')

%% CS_SL0
A=Phi;    %恢复矩阵 
[theta,Spare_L0]=CS_SL0( y,A,0.001,0.9,2,3);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('近似L0范数追踪结果','真实稀疏度')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('SL0恢复信号','原始信号')

%% CS_UALP
A=Phi;    %恢复矩阵 
[theta]=CS_UALP( y,A,0.1);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('UALP恢复信号','原始信号')

%% CS_RSL0
A=Phi;    %恢复矩阵 
[theta,Spare_L0]=CS_RSL0( y,A,0.001,0.9,2,3,0.001);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('近似L0范数追踪结果','真实稀疏度')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RSL0恢复信号','原始信号')


%% CS_IRLS
A=Phi;    %恢复矩阵 
[theta]=CS_IRLS( y,A,0,1e-8,0.1);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('IRLS恢复信号','原始信号')

%% CS_RL1    加权L1最小化，相当于L0最小化
A=Phi;    %恢复矩阵 

[theta]=CS_RL1( y,A,3);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RL1恢复信号','原始信号')

%% CS_IHT    加权L1最小化，相当于L0最小化
A=Phi;    %恢复矩阵 

% [theta]=CS_IHT( y,A,K);
[theta]=CS_IHT( y,A);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('IHT恢复信号','原始信号')


%% CS_SBI    
A=Phi;    %恢复矩阵 

theta=CS_SBIL1(y,A,2000,100,1e4)
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('SBIL1恢复信号','原始信号')

%% ISTA
A=Phi;    %恢复矩阵  
[theta,erro_rnn]=CS_ISTA( y,A,0.00819); %0.00819
figure
plot(erro_rnn,'-*')
legend('ISTA误差')
%% 恢复信号和原始信号比较
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('ISTA恢复信号','原始信号')

%% FISTA
A=Phi;    %恢复矩阵  
[theta,erro_rnn]=CS_ISTA( y,A,0.00819); %0.00819
figure
plot(erro_rnn,'-*')
legend('FISTA误差')
%% 恢复信号和原始信号比较
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('FISTA恢复信号','原始信号')





