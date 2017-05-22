
%%
%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(l1-MAGIC和l1_ls解l1问题)   
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
% Phi=(sqrt(N))*eye(M)*Phi;
Psi=(sqrt(N))*eye(N,N);  %矩阵Psi构造，使得信号x稀疏，由于x本身就是稀疏的，所以这里就是单位采样矩阵

%% 3. 量测信号
y=Phi*x;
%% 4. 重构信号 l1最小化   Using  l1-MAGIC
% A=Phi;   %恢复矩阵,稀疏化矩阵为单位矩阵，因为信号本身就是稀疏的，不需要做任何稀疏变换
A=Phi*Psi';
x0=A'*y;  %最小二乘解估计一个初始值
xh1=l1eq_pd(x0,A,[],y,1e-3);

%% 5. 重构信号l1最小化   Using l1_ls
lambda  = 0.01; % 正则化参数
rel_tol = 1e-3; % 目标相对对偶间隙
quiet=1;   %不输出中间结果
[xh2,status]=l1_ls(A,y,lambda,rel_tol,quiet);

%% 6.恢复信号和原始信号比较
figure
plot(t,Psi*xh1,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-MAGIC恢复信号','原始信号')

figure
plot(t,Psi*xh2,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-ls恢复信号','原始信号')
