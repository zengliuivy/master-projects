clear all; close all;
warning('off')

fig_num =0;
%蒙特卡洛仿真次数
mc=1;

 for h=1:mc

%        tic
% %代码块
% 扫描次数与扫描周期
N = 4;
M = 3;
T = 5; %秒

% 所考虑的正方形仿真区域
Xscope = 10^5;
Yscope = 10^5;

% 目标运动参数
v = 500;     % 500m/s
theta = 0;   % 水平正x轴运动

sigmax=1; %x轴方向的随机加速度
sigmay=0.6; %y轴方向的随机加速度

%%%%%设置噪声%%%%%%%%%%%
S=v/5;%%设置噪声方差
w=sqrt(S);%%设置均方根

% 距离观测标准差与方位角观测标准差
sigma_r = w;%%高斯噪声标准差
sigma_theta = 0.3;

% 所考虑的正方形仿真区域内的杂波平均数
renbuda =90;

% 指定几次扫描的杂波个数，每个周期的数目服从泊松分布，分布的均值由面积大小
% 以及单位面积内杂波数的乘积确定
K = poissrnd(renbuda, 1, N);

% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% 限制关联规则中的最大与最小速度、最大加速度和最小角速度，连续三次扫描的夹角
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
vw=6;%%转弯角速度，单位为度/s，转弯角速度的设置与波门设计有关
thetamax =30;
yama=25;
% thetamax = pi/2;

%量测方程
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
 %P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
 R = [10000  0;0 10000];%%%设置与量测噪声有关
% R = [sigma_r^4  sigma_theta^4;sigma_theta^4 sigma_r^4];
P0=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];


Temp(1).sample=[]; Temp(2).sample=[];Temp(3).sample=[];Temp(4).sample=[];%%%结构体，用于存储杂波过滤后的数据
if ~exist('track','var')
else
for r=1:size(track,2)
track(r).assoi_point=[];track(r).seq=[];
end
end

%% 仿真产生5个目标的航迹(量测数据) %%

radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4行2列
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
%% 每次扫描所得点迹集合sample中的前5个点被设定为目标点 %%
i = 0;
for k = K
    i = i + 1;
     cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycle为结构体   存储杂波
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
    
end

%% 用第一次扫描的点迹建立暂时航迹 %%
for i = 1:size(cycle(1).sample, 1)
    track(i).seq = cycle(1).sample(i,:); %%%存储所有的第一次扫描的结果
%     track(i).shouldadd = [];
    track(k).assoi_point = [];      %存储与航迹关联的点迹的关联值
end

%% 用第二次扫描的点建立可能航迹 %%
for i = 2
    tracknum = size(track,2);      %求得暂态航迹数（即第一次扫描的所有量测数量）
    tracknum_temp = tracknum;
    samplenum = size(cycle(i).sample,1);     %求得第二次扫描所有的量测点迹数
    
    D = zeros(tracknum,samplenum);      %存储暂态航迹与量测的关联值
    %% 计算本次扫描的所有点迹与暂态航迹的关联值 %%
    for j = 1:samplenum
        data = cycle(i).sample(j,:); %%当前第二次扫描的量测值
        for k = 1:tracknum
            if size(track(k).seq,1) > 0 
                data1 = track(k).seq; %%当前第一次扫描的量测值
                D(k,j) = (data(1)-data1(1))^2 + (data(2)-data1(2))^2; %%二者之间的距离平方
                alpha(k,j)=AngX([data(1)-data1(1);data(2)-data1(2)]);%%%目标方向计算
            end                 
        end
    end
    
    for j = 1:samplenum
        flag = 0;
        
        for k = 1:tracknum
            if D(k,j) >= (vmin*T-w)^2 && D(k,j) <= (vmax*T+w)^2 && alpha(k,j)<90%%波门设计,新增目标运动方向限制，沿着x轴正向
                track(k).assoi_point = [track(k).assoi_point;D(k,j) j];
                Temp(1).sample=[Temp(1).sample;track(k).seq(1,1) track(k).seq(1,2)];
                            Temp(2).sample=[Temp(2).sample;cycle(i).sample(j,:)];
                flag = 1;
            end
        end
        
        %% 与暂态航迹未关联的点迹作为新的暂态航迹头
        if flag == 0
            tracknum_temp =tracknum_temp + 1;
            
            track(tracknum_temp).seq = cycle(i).sample(j,:);
            track(tracknum_temp).assoi_point = [];
        end
    end
    
    %% 由关联点迹判别，对暂态航迹进行处理 %%
    for k = 1:tracknum
%         for k = 1:tracknum_temp
        L = size(track(k).assoi_point,1);
        if L == 1 %%如果已经关联
            j = track(k).assoi_point(end,2);%第二次扫描的第几个点迹
            track(k).seq = [track(k).seq;cycle(i).sample(j,:)];%%将关联上的点迹数据都放在一个seq里
        end
        if L > 1%如果有好几个点关联，选择最近的那个点
            min = track(k).assoi_point(1,:);
            for j = 2:L
                if (track(k).assoi_point(j,1) - (v*T)^2) < (min(1) - (v*T)^2)
                    min = track(k).assoi_point(j,:);
                end
            end
            track(k).seq = [track(k).seq;cycle(i).sample(min(2),:)];
        end
        if L == 0 %如果没关联，直接清空
            track(k).seq = [];
        end
    end
    
    %% 整编航迹 %%
    track1 = [];
    track1num = 0;
    for j = 1:tracknum_temp
        if ~isempty(track(j).seq)
            track1num = track1num + 1;
            track1(track1num).seq = track(j).seq;
            track1(track1num).assoi_point = track(j).assoi_point;
        end
    end
    track = track1;
%     clear track1;
end

%% 用第三帧扫描的点迹建立航迹，改用扇形波门 %%
for i = 3
    tracknum = size(track,2);
    samplenum = size(cycle(i).sample,1);
    for j = 1:tracknum
        track(j).assoi_point = [];
        num = size(track(j).seq,1);     %计算每条航迹中的点迹数
        %% 暂态航迹中只有一个点迹时，在第三次量测数据中重复上述操作，依旧使用圆形波门%% 
        if num == 1   
            data = track(j).seq(end,:);
            for k = 1:samplenum
                data1 = cycle(i).sample(k,:);
                d(k) = (data1(1) - data(1)*v*T )^2 + (data1(2) - data(2))^2;%%%%对目标进行线性外推
                alpha(k,j)=AngX([data1(1)-data(1)*v*T ;data1(2)-data(2)]);%%%目标方向计算
                if d(k) >= (vmin*T-w)^2 && d(k) <= (vmax*T+w)^2  && alpha(k,j)<90
                    track(j).assoi_point = [track(j).assoi_point;d(k) k];
                            Temp(1).sample=[Temp(1).sample;track(j).seq(1,1) track(j).seq(1,2)];
                            Temp(3).sample=[Temp(3).sample;cycle(i).sample(k,:)];
                end
            end
            L = size(track(j).assoi_point,1);
            if L == 1
                k = track(j).assoi_point(end,2);
                track(j).seq = [track(j).seq;cycle(i).sample(k,:)];
            end
            if L > 1
               min = track(j).assoi_point(1,:);
               for  k = 2:L
                    if (track(j).assoi_point(k,1) - (v*T)^2) < (min(1) - (v*T)^2)
                        min =track(j).assoi_point(k,:);
                    end
               end
               track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];
            end
            if L == 0
            track(j).seq = [];
            end
        end
        %% 暂态航迹中多于一个点迹时 %%
        if num > 1
            data = track(j).seq(end-1,:);     %航迹中的倒数第二个点迹
            data1 = track(j).seq(end,:);      %航迹中的最后一个点迹
            %%外推
%             X1 = [data1(1) (data1(1)-data(1))/T data1(2) (data1(2)-data(2))/T]';%%该航迹横坐标 x轴速度 纵坐标 Y轴速度
% 
%             Z1 = H*F*X1;%%%外推点坐标
            %卡尔曼滤波改进
            X=[data(1);v;data(2);0];
%         P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
   P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
            x_ =F*X;          %预测当前状态
            P1 = F*P*F';
            % 更新步
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data1(1);data1(2)] - H*x_);
            P= (1-k*H)*P1;
           
            x_ = F*x;            %预测当前状态
            Z1=H*x_;
            P1 = F*P*F';
            S = H*P1*H'+ R;%%新息协方差矩阵
            for k = 1:samplenum 
                data2 = cycle(i).sample(k,:)';%%%%第三圈量测数据
%                 d1(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
%                 a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
%                 b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
%                 c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
%                 alpha = acos((a + b - c)/(2*sqrt(a*b)));
%                 alpha = 180 - alpha * 180 / pi;
                d1(k) =sqrt((data2(1)-data1(1))^2 + (data2(2)-data1(2))^2);%计算量测到b点距离
                BBC=sqrt((Z1(1)-data1(1))^2 + (Z1(2)-data1(2))^2);%计算bc点距离
                bc=[data2(1)-data1(1);data2(2) - data1(2)];%%bc'向量
                bbc=[Z1(1) - data1(1);Z1(2) - data1(2)];%%bc向量
               cc=sqrt((Z1(1)-data2(1))^2+(Z1(2)-data2(2))^2);%cc'之间的距离
%                 %计算机动转弯角
                alpha=abs(AngX(bc)-AngX(bbc));
%                a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
%                b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
%                c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
%                alpha = acos((a + b - c)/(2*sqrt(a*b)));
%                alpha = 180 - alpha * 180 / pi;
%            plot([data(1);data1(1);Z1(1)],[data(2);data1(2);Z1(2)],'o');hold on;
%            plot(radar1(:,1),radar1(:,2));
                %计算扇形区域半径
                ef=(vmax-vmin)*T+((amax-amin)*T^2)/2+w;
                 
                if  abs( d1(k)-BBC)<ef/2 && alpha < vw*T  %%速度波门  
                    track(j).assoi_point = [track(j).assoi_point;cc k];
                      Temp(3).sample=[Temp(3).sample;cycle(i).sample(k,:)];
                end
            end
            L = size(track(j).assoi_point,1);
            if L == 1
                k = track(j).assoi_point(end,2);
                track(j).seq = [track(j).seq;cycle(i).sample(k,:)];%%添加新的量测点到航迹
            end
            if L > 1 %%若不止一个点与其关联
                min = track(j).assoi_point(1,:);
                for k = 2:L
                    if track(j).assoi_point(k,1) < min(1)
                        min = track(j).assoi_point(k,:);
                    end
                end
                track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];%%%选择距离外推点c最近的量测相关联
            end
        end
    end
     
    %% 对航迹进行整编 %% 
    track2 = [];
    track2num = 0;
    for j = 1:tracknum
        if ~isempty(track(j).seq)
            track2num = track2num + 1;
            track2(track2num).seq = track(j).seq;
            track2(track2num).assoi_point = track(j).assoi_point;
        end
    end
    track = track2;
    
end

%% 用第四次扫描的点迹继续判别航迹，采用椭圆波门 %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %第四次扫描的量测数
    tracknum = size(track,2);    %求得此时的航迹数
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %如果航迹在上一帧没有关联点迹，则对航迹进行外推两次
          X =  F*[track(j).seq(end,1) (track(j).seq(end,1)-track(j).seq(end-1,1))/T track(j).seq(end,2) (track(j).seq(end,2)-track(j).seq(end-1,2))/T]'; 
          data1 = H*X;
%           track(j).seq = [track(j).seq;data1'];
           X=[track(j).seq(end,1);v;track(j).seq(end,2);0];
%            P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
%  R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
%     R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
%     R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
           P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            x_ =F*X;          %预测第三帧状态
            P1 = F*P*F';
            % 更新步
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data1(1);data1(2)] - H*x_);
            P= (1-k*H)*P1;
                  
             x_ = F*x; %预测第四帧外推点
            Z1=H*x_;
        P1 = F*P*F';
           S = H*P1*H' +R;%%新息协方差矩阵
          for k = 1:samplenum
              data2 = cycle(i).sample(k,:)';
              d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
              a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
              b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
              c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
              alpha = acos((a + b - c)/(2*sqrt(a*b)));
              alpha = 180 - alpha * 180 / pi;
              if d(k) <yama && alpha < thetamax
                  track(j).assoi_point = [track(j).assoi_point;d(k) k];
                  Temp(4).sample=[Temp(4).sample;cycle(i).sample(k,:)];
              end
          end
          if size(track(j).assoi_point,1) == 1
              k = track(j).assoi_point(end,2);
              track(j).seq = [track(j).seq;cycle(i).sample(k,:)];
          end
          if size(track(j).assoi_point,1) > 1
              min = track(j).assoi_point(1,:);
              for k = 2:size(track(j).assoi_point,1)
                  if track(j).assoi_point(k,1) < min(1)
                      min = track(j).assoi_point(k,:);
                  end
              end
              track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];
          end
          if size(track(j).assoi_point,1) == 0
              track(j).seq = [];
          end
        end
        if L > 0
           track(j).assoi_point = [];
           data = track(j).seq(end-1,:);
           data1 = track(j).seq(end,:);
%            X = [data1(1) (data1(1) - data(1))/T data1(2) (data1(2)-data(2))/T]';
%            P = F*[sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2]*F';
%            X1 = F*X;
%            P1 = F*P*F';%%预测协方差方程
%            Z1 = H*X1;%%预测量测位置矩阵
           
              %卡尔曼滤波改进
            X=[track(j).seq(end-2,1);v;track(j).seq(end-2,2);0];
           P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
     %       P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            x_ =F*X;          %预测第二帧状态
            P1 = F*P*F';
            % 更新步
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data(1);data(2)] - H*x_);
            P= (1-k*H)*P1;
           
            x_ = F*x;            %预测第三帧状态
            P1 = F*P*F';
            %更新
             k = P1*H'/(H*P1*H'+R);
             x = x_ + k*([data1(1);data1(2)] - H*x_);
             P= (1-k*H)*P1;
                  
             x_ = F*x; %预测第四帧外推点
            Z1=H*x_;
            P1 = F*P*F';
           S = H*P1*H' +R;%%新息协方差矩阵
%     X=[data(1);v;data(2);0];
%          P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
% % P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
% %     R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
% %     R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
% %     R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
%             x_ =F*X;          %预测第三帧状态
%             P1 = F*P*F';
%             % 更新步
%             k = P1*H'/(H*P1*H'+R);
%             x = x_ + k*([data1(1);data1(2)] - H*x_);
%             P= (1-k*H)*P1;
%            
%             x_ = F*x;            %预测当前状态
%             Z1=H*x_;
%             S = H*P1*H'+ R;%%新息协方差矩阵
           for k = 1:samplenum
               data2 = cycle(i).sample(k,:)';
               d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
               a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
               b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
               c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
               alpha = acos((a + b - c)/(2*sqrt(a*b)));
               alpha = 180 - alpha * 180 / pi;
%                 plot([data(1);data1(1);Z1(1);data2(1)],[data(2);data1(2);Z1(2);data2(2)],'o');hold on;
%            plot(radar5(:,1),radar5(:,2));
               if d(k) < yama && alpha < thetamax %%参数查询X^2的分布表便可得到
                   track(j).assoi_point = [track(j).assoi_point;d(k) k];
                   Temp(4).sample=[Temp(4).sample;cycle(i).sample(k,:)];
               end
           end
           if size(track(j).assoi_point,1) == 1
               k = track(j).assoi_point(end,2);
               track(j).seq = [track(j).seq;cycle(i).sample(k,:)];
           end
           if size(track(j).assoi_point,1) > 1
              min = track(j).assoi_point(1,:);
              for k = 2:size(track(j).assoi_point,1)
                  if track(j).assoi_point(k,1) < min(1)
                      min = track(j).assoi_point(k,:);
                  end
              end
              track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];
           end
        end
    end
    
%     %% 整编航迹 %%

    tracknew = [];
    tracknewnum = 0;
    for j = 1:tracknum
        if size(track(j).seq,1) > 2
            tracknewnum = tracknewnum + 1;
            tracknew(tracknewnum).seq = track(j).seq;
        end
    end
end

%% 模糊的的hough变换 %%
target=5;%目标数
k=200;%sig分的个数
m=400;%p分的个数
Ln=150000;%雷达的量测距离

Np=1:k;
dNp=pi/k;%参数空间角度间隔
% dNp=180/k;%参数空间角度间隔
angle=(Np-1/2)*dNp;

dMp=2*Ln/m;%参数空间垂距间隔
[~,Seq]=size(tracknew);
for i=1:Seq
A(i).seq=zeros(k,m);%积累矩阵初始化，设为全零矩阵
A1(i).seq=zeros(k,m);%积累矩阵初始化，设为全零矩阵
end
AA=zeros(k,m);
 
%% hough变换 %%
  %计算相应的pn估值
    for i=1: Seq
        L=size(tracknew(i).seq,1);
        for l=1:L
          for j=1:k
            Pn(i).sample(l,j)=tracknew(i).seq(l,1)*cos(angle(j))+tracknew(i).seq(l,2)*sin(angle(j));
          end
        end
    end

    %对积累矩阵投票,看pn是否落在方格里，如果是就加一，设置阈值筛选，主要为了减少计算量
  for l=1: Seq
    for i=1:k
        for j=1:m
%             a=-T+(j-1)*dMp;
%             b=-T+j*dMp;
            a=(j-1)*dMp;
            b=j*dMp;
           for h=1:size(tracknew(l).seq,1)
               if (Pn(l).sample(h,i)>=a && Pn(l).sample(h,i)<b) 
                   A1(l).seq(i,j)=A1(l).seq(i,j)+1;
%                    AA(i,j)=AA(i,j)+1;
               end
           end
        end
    end
  end
   
   %对积累矩阵投票,利用模糊隶属度函数
%   sita_m=deg2rad(90);%%最大误差范围
  sita_m=1;%%最大误差范围
  pm=1400;%%%p向最大误差范围
%   zita_sita=1;
%    zita_p=200;
  index=1;
 for l=1: Seq
      for i=1:k
        for j=1:m
            if A1(l).seq(i,j)>2%%%设置较低的阈值，为了减少计算量
%                 num=1;
% l=3;
                uk=0;
%                  a=-T+(j-1)*dMp;
%                  b=-T+j*dMp;
            a=(j-1/2)*dMp;
            b=angle(i);  
            if i>2 && i<=k-2
                f=i-2;
                e=i+2;
            elseif i<=2
                f=1;
                  e=i+2;
            elseif i>k-2
                 f=i-2;
                  e=k;      
            end
               for h=1:size(tracknew(l).seq,1)
                   for g=f:1:e
                 if abs(Pn(l).sample(h,g)-a)<=2*dMp  %%是否落入方格附近的格子内
                         pii=angle(g);
                         pj=Pn(l).sample(h,g);%将该峰值点作为自变量  
                         for ii=1:(2*sita_m/dNp+1)
                            sita_ijl(ii)=pii-sita_m+(ii-1)*dNp;
                         end
                         for jj=1:(2*pm/dMp+1)
                             p_ijl(jj)=pj-pm+(jj-1)*dMp;
                         end
                          zita_sita=std(sita_ijl,0);%%计算θ 方向上的均方差
                          zita_p=std(p_ijl,0);%%计算 ρ方向上的均方差
                           uk=uk+gaussmf(pii,[ zita_sita b])* gaussmf(pj,[zita_p a]);%%%隶属度函数
%                         A(l).coordinate(i,j)=index;%%建立索引
%                         Co(index).seq(num,:)=[tracknew(l).seq(h,1),h,i,j];%%%用来存储x坐标，便于还原轨迹
%                         num=num+1;
                 end
                   end
               end
                  A(l).seq(i,j)=A(l).seq(i,j)+uk;%%%%累积矩阵的计算
%                   index=index+1;
             end
           end
      end
 end  
%  bar3(A(3).seq);
% for l=1:Seq
%       AA=AA+A(l).seq;
% end
%  bar3(AA)
%     xlabel('ρ');
%     ylabel('θ');
%     zlabel('Cumulative value');
%       max_sum=max(max(A(3).seq));
%      [maxi,maxj]=find(A(3).seq==max_sum);
%         A(3).seq(maxj,maxi)=0;
%      plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
%      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
%  hold on;
%         X=0:1:100000;
%                          P0=(66-1/2)*dMp;
%                   A0=angle(47);  
%         YS=(P0-X*cos(A0))/(sin(A0));
%         plot(X,YS,'k');
%         hold on
 %% 阈值法确认航迹%%
 
 num=1;track_hough=[];
 %%寻找全局最值
 max_int=max(max(A(1).seq));
 jihe(1)=max_int;
 for l=2: Seq
             max_sum=max(max(A(l).seq));
             jihe(l)=max_sum;
             if max_sum>max_int%%
              max_int=max_sum;
             end
 end
 nu=1;
 for l=1: Seq
             if jihe(l)>max_int*0.5%%
                 ji(nu)=jihe(l);
                 nu=nu+1;
             end
 end
 max_mean=mean(ji(:));
 for l=1: Seq
             max_sum=max(max(A(l).seq));
             if max_sum> max_mean*0.68%%设置阈值
%            if max_sum>8%%设置阈值
                [maxi,maxj]=find(A(l).seq==max_sum);
                 P0=(maxj-1/2)*dMp;
                 A0=angle(maxi);
                 x1=tracknew(l).seq(1,1);
                 x2=tracknew(l).seq(end,1);
                 track_hough(num,:)=[P0,A0,x1,x2];
                 num=num+1;
%                  A(l).seq(maxi,maxj)=0;
%                  max_sum=max(max(A(l).seq));
             end
end


%% 画图 %%
fig_num = fig_num+1;
figure(fig_num);
 s = ['*', 's', '+', '.'];
hold on;
for i = 1:N
    plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
end
% fig_num = fig_num+1;
% figure(fig_num);
hold on;
plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
     [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
% % for ii = 1:size(tracknew,2)
% %     data = tracknew(ii).seq;
% %     plot(data(:,1),data(:,2), '-');
% % end 
% for ii = 1:size(track2,2)
%     data = track2(ii).seq;
%     plot(data(:,1),data(:,2), '-');
% end 
% 
 plot(radar1(:,1),radar1(:,2), '-');
 plot(radar2(:,1), radar2(:,2), '-');
 plot(radar3(:,1),radar3(:,2), '-');
 plot(radar4(:,1),radar4(:,2), '-');
 plot( radar5(:,1), radar5(:,2), '-');
xlim([0, Xscope]);
ylim([0, Yscope]);
XLabel = xlabel('X/m');
YLabel = ylabel('Y/m'); 
box on;

% 
fig_num = fig_num+1;
figure(fig_num);
hold on;
  x2=scatter([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)],'k','o');
for i = 1:N
%    c = linspace(10,1,length(Temp(i).sample));
     x1=scatter(Temp(i).sample(:,1), Temp(i).sample(:,2),30,'filled');x1.MarkerFaceColor=[0.67 0.67 1];
 end

xlim([0, Xscope]);
ylim([0, Yscope]);
box on;
legend([x1,x2],'Filtered clutters and targets','True targets');
XLabel = xlabel('X/m');
YLabel = ylabel('Y/m'); 


% % % %% hough变换结果
fig_num = fig_num+1;
figure(fig_num);
s = ['*', 's', '+', '.'];
hold on;
for i = 1:N
    plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
end
 hold on;

plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');

  for  s=1:size(track_hough,1)    
        X=track_hough(s,3):1:track_hough(s,4);
        YS=(track_hough(s,1)-X*cos(track_hough(s,2)))/(sin(track_hough(s,2)));
        plot(X,YS,'k');
        hold on
  end
  XLabel = xlabel('X/m');
  YLabel = ylabel('Y/m'); 
% toc
% disp(['运行时间: ',num2str(toc)]);
 end
% % 



