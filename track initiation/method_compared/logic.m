% %% Track Start  五个目标航迹起始  3/4逻辑

clear all; close all;
warning('off')
fig_num =0;
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
renbuda =120;
    
% 指定几次扫描的杂波个数，每个周期的数目服从泊松分布，分布的均值由面积大小
% 以及单位面积内杂波数的乘积确定
K = poissrnd(renbuda, 1, N);

% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% 限制关联规则中的最大与最小速度、最大加速度和最小角速度，连续三次扫描的夹角
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
thetamax =30;
yama=25;
% thetamax = pi/2;

%量测方程
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
% R = [sigma_r^2  sigma_theta^2;sigma_theta^2 sigma_r^2];%%%设置与量测噪声有关
R = [25000  0;0 25000];
P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
%   P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];

%蒙特卡洛仿真实验次数

mc=1;

for h=1:mc
    tic
%代码块
Temp(1).sample=[]; Temp(2).sample=[];Temp(3).sample=[];Temp(4).sample=[];%%%结构体，用于存储杂波过滤后的数据


%% 仿真产生5个目标的航迹(量测数据) %%
% radar1 = simutrack(55000, 55000, v, theta, 0, 0, sigma_r, sigma_theta, T, N); %4行2列
% radar2 = simutrack(45000, 45000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar3 = simutrack(35000, 35000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar4 = simutrack(45000, 25000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar5 = simutrack(55000, 15000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar1 = simutrack(55000, 55000, v, theta, 0, 0, 0, 0, T, N); %4行2列
% radar2 = simutrack(45000, 45000, v, theta, 0, 0, 0, 0, T, N);
% radar3 = simutrack(35000, 35000, v, theta, 0, 0, 0, 0, T, N);
% radar4 = simutrack(45000, 25000, v, theta, 0, 0, 0, 0, T, N);
% radar5= simutrack(55000, 15000, v, theta, 0, 0, 0, 0, T, N);
% 
 radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4行2列
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
%% 每次扫描所得点迹集合sample中的前5个点被设定为目标点 %%
i = 0;
for k = K
    i = i + 1;
    cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycle为结构体   存储杂波点
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
end

%% 用第一次扫描的点迹建立暂时航迹 %%
for i = 1:size(cycle(1).sample, 1)
    track(i).seq = cycle(1).sample(i,:);
%     track(i).shouldadd = [];
    track(k).assoi_point = [];      %存储与航迹关联的点迹的关联值
end

%% 用第二次扫描的点建立可能航迹 %%
for i = 2
    tracknum = size(track,2);      %求得暂态航迹数
    tracknum_temp = tracknum;
    samplenum = size(cycle(i).sample,1);     %求得第二帧的量测点迹数
    
    D = zeros(tracknum,samplenum);      %存储暂态航迹与量测的关联值
    %% 计算本次扫描的所有点迹与暂态航迹的关联值 %%
    for j = 1:samplenum
        data = cycle(i).sample(j,:);
        for k = 1:tracknum
            if size(track(k).seq,1) > 0
                data1 = track(k).seq;
                D(k,j) = (data(1)-data1(1))^2 + (data(2)-data1(2))^2;
            end                 
        end
    end
    
    for j = 1:samplenum
        flag = 0;
        for k = 1:tracknum
            if D(k,j) >= (vmin*T)^2 && D(k,j) <= (vmax*T)^2
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
        L = size(track(k).assoi_point,1);
        if L == 1
            j = track(k).assoi_point(end,2);
            track(k).seq = [track(k).seq;cycle(i).sample(j,:)];
        end
        if L > 1
            min = track(k).assoi_point(1,:);
            for j = 2:L
                if (track(k).assoi_point(j,1) - (v*T)^2) < (min(1) - (v*T)^2)
                    min = track(k).assoi_point(j,:);
                end
            end
            track(k).seq = [track(k).seq;cycle(i).sample(min(2),:)];
        end
        if L == 0
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

%% 用第三帧扫描的点迹建立可靠航迹 %%
for i = 3
    tracknum = size(track,2);
    samplenum = size(cycle(i).sample,1);
    for j = 1:tracknum
        track(j).assoi_point = [];
        num = size(track(j).seq,1);     %计算每条航迹中的点迹数
        %% 暂态航迹中只有一个点迹时 %% 
        if num == 1   
            data = track(j).seq(end,:);
            for k = 1:samplenum
                data1 = cycle(i).sample(k,:);
                d(k) = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
                if d(k) >= (vmin*T)^2 && d(k) <= (vmax*T)^2
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
            X = [data1(1) (data1(1)-data(1))/T data1(2) (data1(2)-data(2))/T]';
%             P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            X1 = F*X;
            P1 = F*P*F';
            Z1 = H*X1;
            S = H*P1*H' + R;
            for k = 1:samplenum 
                data2 = cycle(i).sample(k,:)';
%                 S =[10500,0;0,1.0000036];
                d1(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
                a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
                b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
                c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
                alpha = acos((a + b - c)/(2*sqrt(a*b)));
                alpha = 180 - alpha * 180 / pi;
                if d1(k) < yama && alpha < thetamax
                    track(j).assoi_point = [track(j).assoi_point;d1(k) k];
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
                for k = 2:L
                    if track(j).assoi_point(k,1) < min(1)
                        min = track(j).assoi_point(k,:);
                    end
                end
                track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];
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
    
%     %% 对于可靠航迹中的点迹有等于三个的，作为成功起始的新航迹 %%
%     tracknew = [];
%     tracknewnum = 0;
%     for j = 1:track2num
%         if size(track(j).seq,1) == 3
%             tracknewnum = tracknewnum + 1;
%             
%             tracknew(tracknewnum).seq = track(j).seq;
%         end
%     end
end

%% 用第四次扫描的点迹继续判别航迹 %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %第四次扫描的量测数
    tracknum = size(track,2);    %求得此时的航迹数
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %如果航迹在上一帧没有关联点迹，则对航迹进行直线外推
          X =  F*[track(j).seq(end,1) (track(j).seq(end,1)-track(j).seq(end-1,1))/T track(j).seq(end,2) (track(j).seq(end,2)-track(j).seq(end-1,2))/T]'; 
          P0 = F*P*F';
          X1 = F*X;
          P1 = F*P0*F';
          Z1 = H*X1;
          S = H*P1*H' + R;
          data1 = [data(1)*v*T data(2)];
%           track(j).seq = [track(j).seq;data1];
          for k = 1:samplenum
              data2 = cycle(i).sample(k,:)';
              d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
              a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
              b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
              c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
              alpha = acos((a + b - c)/(2*sqrt(a*b)));
              alpha = 180 - alpha * 180 / pi;
              if d(k) < yama && alpha < thetamax
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
           X = [data1(1) (data1(1) - data(1))/T data1(2) (data1(2)-data(2))/T]';
           P0 = F*P*F';
           X1 = F*X;
           P1 = F*P0*F';
           Z1 = H*X1;
           S = H*P1*H' + R;
           for k = 1:samplenum
               data2 = cycle(i).sample(k,:)';
               d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
               a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
               b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
               c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
               alpha = acos((a + b - c)/(2*sqrt(a*b)));
               alpha = 180 - alpha * 180 / pi;
               if d(k) < yama && alpha < thetamax
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
%     track3 = [];
%     track3num = 0;
%     for j = 1:tracknum
%         if ~isempty(track(j).seq)
%             track3num = track3num + 1;
%             
%             track3(track3num).seq = track(j).seq;
%         end
%     end
    %%  对于可靠航迹中的点迹有等于三个的，作为成功起始的新航迹 %%   
%     for j = 1:track3num
%         if size(track3(j).seq,1) == 3
%             tracknewnum = tracknewnum + 1;
%             
%             tracknew(tracknewnum).seq = track3(j).seq;
%         end
%     end
    tracknew = [];
    tracknewnum = 0;
    for j = 1:tracknum
        if size(track(j).seq,1) > 2
            tracknewnum = tracknewnum + 1;
            tracknew(tracknewnum).seq = track(j).seq;
        end
    end
end

fig_num = fig_num+1;
 figure(fig_num);
% s = ['*', 's', '+', '.'];
hold on;
% for i = 1:N
%     plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
% end

plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
     [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
for ii = 1:size(tracknew,2)
    data = tracknew(ii).seq;
    plot(data(:,1),data(:,2), '-');
end 
xlim([0, Xscope]);
ylim([0, Yscope]);
box on;
% toc
% disp(['运行时间: ',num2str(toc)]);

% fig_num = fig_num+1;
% figure(fig_num);
% hold on;
% for i = 1:N
% %    c = linspace(10,1,length(Temp(i).sample));
% if isempty(Temp(i).sample)
%     Temp(i).sample=[0 0];
% end
%      x1=scatter(Temp(i).sample(:,1), Temp(i).sample(:,2),30,'filled');x1.MarkerFaceColor=[0.67 0.67 1];
% end
%  x2=scatter([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
%       [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)],'k','o');
% % for s=1:targets     
% %         X=0:1:1000;
% %         YS=(P0(s)-X*cos(A0(s)))/(sin(A0(s)));
% %         plot(X,YS,'k');
% %         hold on
% % end
% xlim([0, Xscope]);
% ylim([0, Yscope]);
% box on;
% legend([x1,x2],'Filtered clutters and targets','True targets');
end
