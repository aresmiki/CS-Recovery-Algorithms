clc;clear all;close all;
%%
%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(l1-MAGIC和l1_ls解l1问题)   
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--西南交通大学牵引动力国家重点实验室 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月26日
%---------------------------------------------------------------------------------%
%% 1. 构造信号
fs=100;     %采样频率
N=100;    %信号长度
t=0:1/fs:(N-1)/fs; 
x2=cos(2*pi*50*t);  %构造信号
%% 2. 离散余弦变换，并将很小值设置为0，确保稀疏度，并重构回信号
% C=gen_dct(N);
C=dctmtx(N);     %离散余弦变换矩阵
cx=C*x2';
cx(find(abs(cx)<0.5))=0;   %将较小分量置零，虽然影响原始信号，但确保了稀疏度
% figure
% plot([x2',C'*cx])
x2=C'*cx;    %重构回信号，该信号的离散余弦必定稀疏
x2=x2';
%% 3. 测量信号   
% 用44测量信号的数据恢复100个点的数据，按照奈奎斯特Nyquist采样定理，1s需采样100个点才能恢复原始信号，
% 但是CS理论只需要44个点的数据就能恢复，这完全突破了Nyquist采样定理的限制。

K=length(find(abs(cx)>0.5));   %信号稀疏度,查看离散余弦变换的图
M=2*ceil(K*log(N/K)); %K=9是，该值为22，测量数,测量矩阵压缩程度，经验公式
randn('state',4)
Phi=randn(M,N);  %  测量矩阵(高斯分布白噪声)
Phi=orth(Phi')';    %正交化
y=Phi*x2.';     %  获得线性测量 ---只有44个点，

%% 4. l1问题最小化 l1-Magic工具箱 
% l1eq_pd方法解的问题是如下优化问题
% min_x ||x||_1  s.t.  Ax = b
%
A=Phi*C';  
% x0=A'*y;   %最小二乘作为l1最小化的初始值估计
% 用l1-MAGIC的MATLAB工具箱解l1最小化问题
% xh1=l1eq_pd(x0,A,[],y,1e-3);
xh1=l1eq_pd(zeros(N,1),A,[],y,1e-3);  %可以不给初始的估计
%%  l1问题最小化  l1_ls工具箱
% l1_ls解的优化问题是
% minimize ||A*x-y||^2 + lambda*sum|x_i|,
%
lambda  = 0.01; % 正则化参数
rel_tol = 1e-3; % 目标相对对偶间隙
quiet=1;   %不输出中间结果
[xh2,status]=l1_ls(A,y,lambda,rel_tol,quiet);
% At=A';
% [xh2,status]=l1_ls(A,At,M,N,y,lambda,rel_tol,quiet);
%% 5.恢复信号和原始信号比较
figure
plot(t,C'*xh1,'k.-',t,x2,'r-')
xlim([0,t(end)])
legend('l1-MAGIC恢复信号','原始信号')

figure
plot(t,C'*xh2,'k.-',t,x2,'r-')
xlim([0,t(end)])
legend('l1-ls恢复信号','原始信号')

%% 6. ISTA方法解l1最小化问题
A=Phi*C';  
[theta,erro_rnn]=CS_ISTA( y,A,0.05);
figure
plot(erro_rnn,'-*')
legend('ISTA误差')
%% 7.恢复信号和原始信号比较
figure
plot(t,C'*theta,'k.-',t,x2,'r-')
xlim([0,t(end)])
legend('ISTA恢复信号','原始信号')

%% 8. FISTA方法解l1最小化问题
A=Phi*C';  
[theta,erro_rnn]=CS_FISTA( y,A,0.05);
figure
plot(erro_rnn,'-*')
legend('FISTA误差')
%% 9.恢复信号和原始信号比较
figure
plot(t,C'*theta,'k.-',t,x2,'r-')
xlim([0,t(end)])
legend('FISTA恢复信号','原始信号')
%% 10.解L0问题，用近似方法，也就是Lp方法  SL0算法
% A*x=y, using  Smoothed L0
% minimize ||x||_0
% subject to A*x=y
A=Phi*C';  
[theta]=CS_SL0( y,A,0.001,0.9,1,4);
figure
plot(t,C'*theta,'k.-',t,x2,'r-')
xlim([0,t(end)])
legend('SL0恢复信号','原始信号')
