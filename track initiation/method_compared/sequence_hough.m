%序列Hough变换
clear all; close all;
warning('off')

fig_num =0;
% 扫描次数与扫描周期
N = 12;
M = 3;
T = 5; %秒

% 所考虑的正方形仿真区域
Xscope = 10^5;
Yscope = 10^5;

% 目标运动参数
v = 500;     % 500m/s
theta =0;   % 水平正x轴运动

sigmax=1; %x轴方向的随机加速度
sigmay=0.6; %y轴方向的随机加速度

%%%%%设置噪声%%%%%%%%%%%
S=v/5;%%设置噪声方差
w=sqrt(S);%%设置均方根

% 距离观测标准差与方位角观测标准差
sigma_r = w;%%高斯噪声标准差
sigma_theta = 0.3;

% 所考虑的正方形仿真区域内的杂波平均数
renbuda =30;

%蒙特卡洛仿真实验次数
mc=1;

for h=1:mc
% 指定几次扫描的杂波个数，每个周期的数目服从泊松分布，分布的均值由面积大小
% 以及单位面积内杂波数的乘积确定
K = poissrnd(renbuda, 1, N);

cycle1=[];
%% 仿真产生5个目标的航迹(量测数据) %%

 radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4行2列
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
% tic
% %代码块
%% 每次扫描所得点迹集合sample中的前5个点被设定为目标点 %%
i = 0;
for k = K
    i = i + 1;
     cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycle为结构体   存储杂波
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
end
for l=1:N
 cycle1=[cycle1;cycle(l).sample];
end

target=5;%目标数
targets=6;%序列数
k=200;%sig分的个数
m=400;%p分的个数
L=150000;%雷达量测距离
Pd=1;%检测概率

Np=1:k;
dNp=pi/k;%参数空间角度间隔
angle=(Np-1/2)*dNp;

dMp=2*L/m;%参数空间垂距间隔

A=zeros(k,m);%积累矩阵
 coordinate=zeros(k,m);
% for monte=1:Monte_Carlo
ii=1;max_sum=10;
% while max_sum>targets
for s=1:targets
  %所有量测的 Hough变换
  P=[];
    for i=1:size(cycle1,1)
        for j=1:k
            P(i,j)=cycle1(i,1)*cos(angle(j))+cycle1(i,2)*sin(angle(j));
        end
    end

 
    %对积累矩阵投票
     index=1;
    for i=1:k
        for j=1:m
            num=1;
%             a=L+(j-1)*dMp;
%             b=L+j*dMp;
            a=(j-1)*dMp;
            b=j*dMp;
           for h=1:size(P,1)
               if (P(h,i)>=a && P(h,i)<b) 
                   A(i,j)=A(i,j)+1;
                   coordinate(i,j)=index;%%建立索引
                   Co(index).seq(num,:)=[cycle1(h,1),h,i,j];%%%用来存储x坐标，便于还原轨迹
                   num=num+1;
                   
               end
           end
           index=index+1;
        end
    end
   
         max_init=A(1,1);
        maxi=1;maxj=1;
        for i=1:k
            for j=1:m
                if A(i,j)>=max_init
                    max_init=A(i,j);
                    maxi=i;
                    maxj=j;
                end
            end
        end    
     P0(s)=(maxj-1/2)*dMp;
     A0(s)=angle(maxi);
     hough_index=coordinate(maxi,maxj);
     x1(s)=min(Co(hough_index).seq(:,1));
     x2(s)=max(Co(hough_index).seq(:,1));
     to_ij(s,:)=[maxi,maxj];
        for i=1:size(Co(hough_index).seq,1)
           cycle1(Co(hough_index).seq(i,2),:)=[0,0];
        end
     cycle1(all(cycle1==0,2),:)=[];%%%删除最大峰值所对应量测
 end
    
    %绘图
 fig_num = fig_num+1;
figure(fig_num);
% s = ['*', 's', '+', '.','d','x','>'];
hold on;
for i = 1:N
    plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2),'*');
end
plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
     [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
 hold on;
  for  s=1:size(P0,2)    
        X=x1(s):1:x2(s);
%  X=1:1:100000;
        YS=(P0(s)-X*cos(A0(s)))/(sin(A0(s)));
        plot(X,YS,'k');
        hold on
  end
%     fprintf('总共建立的航迹数量为 %f\n',s); 
%   toc
% disp(['运行时间: ',num2str(toc)]);  
end

        