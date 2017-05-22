%% 含有噪声
% minimize ||x||_1
% subject to: (||Ax-y||_2)^2<=eps;
% minimize :  (||Ax-y||_2)^2+lambda*||x||_1
% y传输中可能含噪 y=y+w
%
%%
clc;clear all;
%% 1.构造一个两个谐波信号
lam=0.37;
itrs=400;
m=380;
sig=0.5;
n=1024;
dt=1/2000;
T=1023*dt;
t=0:dt:T;
t=t(:);
x=sin(697*pi*t)+sin(1975*pi*t);
Dn=dctmtx(n);

%% 2.构造测量矩阵 
rand('state',15);
q=randperm(n);
q=q(:);
y=x(q(1:m));
randn('state',7)
w=sig*randn(m,1);  %产生噪声
yn=y+w;  %压缩矩阵有噪声
Psi1=Dn';
%% 3. 重构信号  ISTA
A=Psi1(q(1:m),:);
% [xh,err]=CS_ISTA(yn,A,lam,itrs);  %ISTA
[xh,err]=CS_ISTA(yn,A,lam);  %ISTA
xx=Psi1*xh;
figure
plot(err,'*-')
legend('ISTA误差')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','ISTA重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','ISTA重构信号')

%% 3. 重构信号  FISTA
A=Psi1(q(1:m),:);
% [xh,err]=CS_FISTA(yn,A,lam,itrs);  %FISTA
[xh,errr]=CS_FISTA(yn,A,lam);  %FISTA
xx=Psi1*xh;
figure
plot(errr,'*-')
legend('FISTA误差')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','FISTA重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','FISTA重构信号')

%% 4. 重构信号  OMP
A=Psi1(q(1:m),:);
[xh,errr]=CS_OMP(yn,A,100);  %OMP
xx=Psi1*xh';
figure
plot(errr,'*-')
legend('OMP误差')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','OMP重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','OMP重构信号')

%% 4. 重构信号  CoSaMP
A=Psi1(q(1:m),:);
[xh,errr]=CS_CoSaMP(yn,A,100);  %CoSaMP
xx=Psi1*xh;
figure
plot(errr,'*-')
legend('CoSaMP误差')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','CoSaMP重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','CoSaMP重构信号')

%% 5. l1问题最小化 l1-Magic工具箱  这个稍差，模型为min_x ||x||_1  s.t.  Ax = b，本身不适合该模型
A=Psi1(q(1:m),:); 
% x0=A'*y;   %最小二乘作为l1最小化的初始值估计
% 用l1-MAGIC的MATLAB工具箱解l1最小化问题
% xh1=l1eq_pd(x0,A,[],y,1e-3);
xh=l1eq_pd(zeros(n,1),A,[],y,1e-3);  %可以不给初始的估计
xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','l1-Magic重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','l1-Magic重构信号')

%% 5. l1问题最小化 l1-ls工具箱   这个模型minimize ||A*x-y||^2 + lambda*sum|x_i|,非常符合该问题
A=Psi1(q(1:m),:); 
lambda  = 0.01; % 正则化参数
rel_tol = 1e-3; % 目标相对对偶间隙
quiet=1;   %不输出中间结果
[xh,status]=l1_ls(A,y,lambda,rel_tol,quiet);
xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','l1-ls重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','l1-ls重构信号')



%% 6. SL0   这个模型minimize （A*x-y） + lambda*sum|x_i|
A=Psi1(q(1:m),:); 

[xh,Spare_L0]=CS_SL0( y,A,0.001,0.9,2,3);

xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-稀疏信号','SL0重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','SL0重构信号')

%% 6. SL0   这个模型minimize ||A*x-y||_2 + lambda*sum|x_i|
A=Psi1(q(1:m),:); 

[xh,Spare_L0]=CS_RSL0( y,A,0.001,0.9,2,3,0.01);

xx2=Psi1*xh;
figure
plot(t,x,'b',t,xx2,'r');
legend('DCT-稀疏信号','RSL0重构信号')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx2(50:100),'r','linewidth',1.5)
legend('DCT-稀疏信号','RSL0重构信号')

figure
plot(t,x,'b',t,xx,'r',t,xx2,'k');
legend('DCT-稀疏信号','SL0重构信号','RSL0重构信号')

figure
bar([norm(x-xx),norm(x-xx2)],'r')



