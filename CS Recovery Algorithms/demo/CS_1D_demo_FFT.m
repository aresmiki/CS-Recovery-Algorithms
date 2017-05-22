%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(正交匹配追踪法Orthogonal Matching Pursuit)   
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人： 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月26日  西南交通大学牵引动力国家重点实验室
%                                        SWJTU  TPL
%---------------------------------------------------------------------------------%
clc
clear all
close all
%% 1. 生成原始信号
fs=400;     %采样频率
f1=25;         %第一个信号频率
f2=50;      %第二个信号频率
f3=100;     %第三个信号频率
f4=200;    %第四个信号频率
N=1024;    %信号长度
t=0:1/fs:(N-1)/fs;   
% x=0.3*cos(2*pi*f1*t)+0.6*cos(2*pi*f2*t)+0.1*cos(2*pi*f3*t)+0.9*cos(2*pi*f4*t);  %构造信号
x=cos(2*pi*f1*t)+cos(2*pi*f2*t)+cos(2*pi*f3*t)+cos(2*pi*f4*t);  %构造信号

%% 1.1查看时域和傅里叶谱
% fx=abs(fftshift(fft(x)))*2/N;
% fsf=(fs/N)*((1:N)-N/2-1);
% figure
% plot(fsf,fx)
%% 可以用下面代码直接将很小的值设置为0
% fft_x=fft(x);
% fft_x(find(abs(fft_x)*2/N<0.1))=0;
% figure
% plot(fsf,fx,fsf,fftshift(fft_x*2/N),'--')
% xx=real(ifft(fft_x));
% figure
% plot(t,x,t,xx,'--')
% x=xx;
%% 2. 时域信号压缩传感，获取测量值
K=8;   %信号稀疏度，傅里叶谱中看出来
M=ceil(K*log(N/K));  %测量数,测量矩阵压缩程度，经验公式
randn('state',2)
Phi=randn(M,N);  %  测量矩阵(高斯分布白噪声)
Phi=orth(Phi')';    %正交化
y=Phi*x';     %  获得线性测量 

%% 3. L_1范数最优化重构信号（有测量值y重构x）
Psi=fft(eye(N,N))/sqrt(N);    %  傅里叶正变换矩阵,认为信号x在傅里叶字典上稀疏：theta=Psi*x.  则：x=Psi'*theta.
% 最小化问题 minimize:       ||theta||_0;
%                  subject to:     y=Phi*Psi'*theta;     ==   令 A=Phi*Psi'.   
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.

%%  4. 正交匹配追踪重构信号
[ fft_y,erro_rn ] = CS_OMP( y,A,2*K );
figure
plot(erro_rn,'-*')
legend('OMP正交匹配追踪误差')
r_x=real(Psi'*fft_y');                         %  做逆傅里叶变换重构得到时域信号

%% 5. 恢复信号和原始信号对比

figure;
hold on;
plot(t,r_x,'k.-')                                 %  重建信号
plot(t,x,'r')                                       %  原始信号
xlim([0,t(end)])
legend('OMP恢复信号','原始信号')


%% 6. CoSaMP 采用压缩采样匹配的方法解l1最小化问题
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.
[theta,erro_rnn]=CS_CoSaMP( y,A,K );
figure
plot(erro_rnn,'-*')
legend('CoSaMP压缩采样追踪误差')
%% 重构并对比
cor_x=real(Psi'*theta);                         %  做逆傅里叶变换重构得到时域信号
figure;
hold on;
plot(t,cor_x,'k.-')                                 %  重建信号
plot(t,x,'r')                                       %  原始信号
xlim([0,t(end)])
legend('CoSaMP恢复信号','原始信号')


%% ISTA和FISTA两种算法解的模型为：
% minimize   ||x||_1
% subject to:||Ax-y||_2<=eps
% second-order cone program(SOCP)，eps是噪声强度  
% 测量矩阵y含有噪声，要估计x就是解上面的凸优化问题
% 无噪声情况是下面的l1问题
% minimize ||x||_1
% subject to: Ax=y;
%% 7. ISTA方法解l1最小化问题---有最优化方法做
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.
[theta,erro_rnn]=CS_ISTA( y,A,0.1,2000);
figure
plot(erro_rnn,'-*')
legend('ISTA误差')
%% 8.恢复信号和原始信号比较
figure
plot(t,real(Psi'*theta),'k.-',t,x,'r-')
xlim([0,t(end)])
legend('ISTA恢复信号','原始信号')

%% 9. ISTA方法解l1最小化问题   有最优化方法做
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.
[theta,erro_rnn]=CS_FISTA( y,A,0.1,2000);
figure
plot(erro_rnn,'-*')
legend('FISTA误差')
%% 10.恢复信号和原始信号比较
figure
plot(t,real(Psi'*theta),'k.-',t,x,'r-')
xlim([0,t(end)])
legend('FISTA恢复信号','原始信号')