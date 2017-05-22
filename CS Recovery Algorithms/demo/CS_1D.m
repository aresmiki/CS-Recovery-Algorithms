%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(正交匹配追踪法Orthogonal Matching Pursuit)   
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--西南交通大学牵引动力国家重点实验室 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月26日
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
Phi=randn(M,N);  %  测量矩阵(高斯分布白噪声)
Phi=orth(Phi')';    %正交化
y=Phi*x';     %  获得线性测量 

%% 3. L_1范数最优化重构信号（有测量值y重构x）
Psi=fft(eye(N,N))/sqrt(N);    %  傅里叶正变换矩阵,认为信号x在傅里叶字典上稀疏：theta=Psi*x.  则：x=Psi'*theta.
% 最小化问题 minimize:       ||theta||_0;
%                  subject to:     y=Phi*Psi'*theta;     ==   令 A=Phi*Psi'.   
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.

%%  4. 正交匹配追踪重构信号
m=2*K;                  %  算法迭代次数(m>=K)
fft_y=zeros(1,N);   %  待重构的谱域(变换域)向量    
Base_t=[];              %  记录基向量的矩阵
r_n=y;                  %  残差值
figure
for times=1:m;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(A(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，val值，pos位置
    Base_t=[Base_t,A(:,pos)];                       %  矩阵扩充，记录最大投影的基向量
    A(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Base_t'*Base_t)^(-1)*Base_t'*y;   %  最小二乘,使残差最小
    r_n=y-Base_t*aug_y;                            %  残差
    erro_rn(times)=norm(r_n,2);
    plot(erro_rn,'r-*')                                        %迭代误差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    if erro_rn(times)<1e-6 %
            break; %跳出for循环
    end
end
legend('OMP正交匹配追踪误差')
fft_y(pos_array)=aug_y;                           %  重构的谱域向量
r_x=real(Psi'*fft_y');                         %  做逆傅里叶变换重构得到时域信号

%% 5. 恢复信号和原始信号对比

figure;
hold on;
plot(t,r_x,'k.-')                                 %  重建信号
plot(t,x,'r')                                       %  原始信号
xlim([0,t(end)])
legend('OMP恢复信号','原始信号')
norm(r_x.'-x)/norm(x)                      %  重构误差

%% 6. CoSaMP 采用压缩采样匹配的方法解l1最小化问题
A=Phi*Psi';                         %  恢复矩阵(测量矩阵*正交反变换矩阵);   x=Psi'*theta.
%function xh=CS_CoSaMP( y,A,K );
    [m,n] = size(y);
    if m<n
        y = y'; %y should be a column vector
    end
    [M,N] = size(A); %传感矩阵A为M*N矩阵
    theta = zeros(N,1); %用来存储恢复的theta(列向量)
    pos_num = []; %用来迭代过程中存储A被选择的列序号
    res = y; %初始化残差(residual)为y
    figure
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
        plot(erro_res,'r-*')                                        %迭代误差
        if norm(res)<1e-6 %Repeat the steps until r=0
            break; %跳出for循环
        end
    end
    theta(pos_num)=theta_ls; %恢复出的theta
%  end
legend('CoSaMP压缩采样追踪误差')
%% 重构并对比
cor_x=real(Psi'*theta);                         %  做逆傅里叶变换重构得到时域信号
figure;
hold on;
plot(t,cor_x,'k.-')                                 %  重建信号
plot(t,x,'r')                                       %  原始信号
xlim([0,t(end)])
legend('CoSaMP恢复信号','原始信号')
norm(r_x.'-x)/norm(x)                      %  重构误差



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
A=Phi*C';  
% x0=A'*y;   %最小二乘作为l1最小化的初始值估计
% 用l1-MAGIC的MATLAB工具箱解l1最小化问题
% xh1=l1eq_pd(x0,A,[],y,1e-3);
xh1=l1eq_pd(zeros(N,1),A,[],y,1e-3);  %可以不给初始的估计
%%  l1问题最小化  l1_ls工具箱
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


%%
%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(l1-MAGIC和l1_ls解l1问题)   
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--西南交通大学牵引动力国家重点实验室 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月27日
%---------------------------------------------------------------------------------%
clc;clear all;
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


%%
%----------------------------------------------------------------------------------%
%  1-D信号压缩传感的实现(l1-MAGIC和l1_ls解l1问题)     信号本身就是稀疏的，
%  不需要稀疏矩阵，恢复矩阵A是单位正交矩阵，用OMP方法求解l1问题.
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--西南交通大学牵引动力国家重点实验室 何刘  Email: aresmiki@163.com
%  编程时间：2017年04月27日
%---------------------------------------------------------------------------------%
clc;clear all;%close all
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
%% 4. 重构信号 l1最小化   Using  l1-MAGIC
x0=A'*y;  %最小二乘解估计一个初始值
xh1=l1eq_pd(x0,A,[],y,1e-3);

%% 5. 重构信号l1最小化   Using l1_ls
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
m=2*K;                  %  算法迭代次数(m>=K)
xh=zeros(1,N);   %  重构向量   
Base_t=[];              %  记录基向量的矩阵
r_n=y;                  %  残差值
A=Phi;    %恢复矩阵
figure
for times=1:m;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(A(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，val值，pos位置
    Base_t=[Base_t,A(:,pos)];                       %  矩阵扩充，记录最大投影的基向量
    A(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Base_t'*Base_t)^(-1)*Base_t'*y;   %  最小二乘,使残差最小
    r_n=y-Base_t*aug_y;                            %  残差
    erro_rn(times)=norm(r_n,2);
    plot(erro_rn,'r-*')                                        %迭代误差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
   if erro_rn(times)<1e-6 %
            break; %跳出for循环
    end
end
legend('OMP正交匹配追踪误差')
xh(pos_array)=aug_y;
figure
plot(t,xh,'ko',t,x,'r.')
xlim([0,t(end)])
legend('OMP恢复信号','原始信号')

%% 8. 用压缩采样匹配追踪(CoSaMP)的方法计算l1最优化问题
%Needell D, Tropp J A. CoSaMP: Iterative signal recovery from incomplete 
% and inaccurate samples [J]. Applied & Computational Harmonic Analysis, 2008, 26(3):301-321.
% 一次性选择2*K个较大的基，每次循环不断删除和补进一定数目的基，最快达到给定精度

m=2*K;                  %  算法迭代次数(m>=K)
theta=zeros(1,N);   %  重构向量   
Base_t=[];              %  记录基向量的矩阵
A=Phi;    %恢复矩阵

%% CoSaMP
%function xh=CS_CoSaMP( y,A,K );
    [m,n] = size(y);
    if m<n
        y = y'; %y should be a column vector
    end
    [M,N] = size(A); %传感矩阵A为M*N矩阵
    theta = zeros(N,1); %用来存储恢复的theta(列向量)
    pos_num = []; %用来迭代过程中存储A被选择的列序号
    res = y; %初始化残差(residual)为y
    figure
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
        plot(erro_res,'r-*')                                        %迭代误差
        if norm(res)<1e-6 %Repeat the steps until r=0
            break; %跳出for循环
        end
    end
    theta(pos_num)=theta_ls; %恢复出的theta
%  end
legend('CoSaMP压缩采样追踪误差')
 %%    
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('CoSaMP恢复信号','原始信号')





%%

clc;clear all;

%% A Proximal-Gradient Algorithm Method

%  minimize ||A*x-y||^2 + lambda*||x||_1     (1)
%  ||A*x-y||^2是二范数平方
%  f(x)= ||A*x-y||^2=x'A'Ax-2x'A'y+y'y
%  一阶导数为：▽f(x)= 2A'Ax-2A'y=2A'(Ax-y)

% (1)的最小化问题转换为 x_k=argmin{(1/(2*t_k))||x-(x_k-t_k▽f(x_k-1))||^2+lambda||x||_1}
% 等价为x_k=argmin{(1/(2*t_k))||x-c_k||^2+lambda||x||_1}
% 其中c_k=x_k-t_k▽f(x_k-1)=x_k-1-2t_kA'(Ax_k-1-y)

% t_k为求解步长
%%




