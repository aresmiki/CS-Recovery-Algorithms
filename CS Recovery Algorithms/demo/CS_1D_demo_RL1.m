%% CS_RL1_test
%------------------------------------------------------------------------------------------%
% The weighted L1 minimization can be viewed as a relaxation of a weighted L0 minimization problem
%  minimize W||x||_0
%  subject to Ax=y
%  编程人： 何刘                                    Email: aresmiki@163.com
%  编程时间：2017年05月3日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%------------------------------------------------------------------------------------------%
clc;clear all;close all
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
% M=2*ceil(K*log(N/K));
M=152;
Phi=randn(M,N);  %高斯矩阵作为感知矩阵
Phi=orth(Phi')';  %正交化

%% 3. 量测信号
y=Phi*x;
A=Phi;%恢复矩阵,稀疏化矩阵为单位矩阵，因为信号本身就是稀疏的，不需要做任何稀疏变换
%% CS_RL1    加权L1最小化，相当于L0最小化
[theta]=CS_RL1( y,A,1);
figure
subplot(3,1,1)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off
[theta]=CS_RL1( y,A,2);
subplot(3,1,2)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off

[theta]=CS_RL1( y,A,5);
subplot(3,1,3)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RL1恢复信号','原始信号')