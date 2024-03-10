clear all; close all;
warning('off')

% 扫描次数与扫描周期
N = 4;
M = 3;
T = 5; %秒

% 所考虑的正方形仿真区域
Xscope = 10^3;
Yscope = 10^3;

% 目标运动参数
v = 5;     % 5m/s
theta = 0;   % 水平正x轴运动

sigmax=0.5; %x轴方向的随机加速度
sigmay=0.06; %y轴方向的随机加速度

%%%%%设置噪声%%%%%%%%%%%
S=v/5;%%设置噪声方差
w=sqrt(S);%%设置均方根

% 距离观测标准差与方位角观测标准差
sigma_r =w;%%高斯噪声标准差
sigma_theta = 0.3;

% 所考虑的正方形仿真区域内的杂波平均数
renbuda = 30;
    
% 指定几次扫描的杂波个数，每个周期的数目服从泊松分布，分布的均值由面积大小
% 以及单位面积内杂波数的乘积确定
K = poissrnd(renbuda, 1, N);


% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% 限制关联规则中的最大与最小速度、最大加速度和最小角速度，连续三次扫描的夹角
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
vw=6;%%转弯角速度，单位为度/s
thetamax = 30;
% thetamax = pi/2;

%量测方程
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
R = [160 0;0 250];


%% 仿真产生5个目标的航迹(量测数据) %%
radar1 = simutrack(550, 550, v, theta, 0, 0, sigma_r, sigma_theta, T, N); %4行2列
radar2 = simutrack(450, 450, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar3 = simutrack(350, 350, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar4 = simutrack(450, 250, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar5 = simutrack(550, 150, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar1 = simutrack(550, 550, v, theta, 0, 0, 0, 0, T, N); %4行2列
% radar2 = simutrack(450, 450, v, theta, 0, 0, 0, 0, T, N);
% radar3 = simutrack(350, 350, v, theta, 0, 0, 0, 0, T, N);
% radar4 = simutrack(450, 250, v, theta, 0, 0, 0, 0, T, N);
% radar5 = simutrack(550, 150, v, theta, 0, 0, 0, 0, T, N);

%  radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4行2列
% radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
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
                d(k) = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
                alpha(k,j)=AngX([data1(1)-data(1);data1(2)-data(2)]);%%%目标方向计算
                if d(k) >= (vmin*T-w)^2 && d(k) <= (vmax*T+w)^2  && alpha(k,j)<90
                    track(j).assoi_point = [track(j).assoi_point;d(k) k];
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
            %外推计算
            X = [data1(1) (data1(1)-data(1))/T data1(2) (data1(2)-data(2))/T]';%%该航迹横坐标 x轴速度 纵坐标 Y轴速度
            P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            X1 = F*X;
            P1 = F*P*F';
            Z1 = H*X1;%%%外推点坐标


            S = H*P1*H' + R;%%新息协方差矩阵
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
                %计算机动转弯角
                alpha=abs(AngX(bc)-AngX(bbc))*2;
            
                %计算扇形区域半径
                ef=(vmax-vmin)*T+((amax-amin)*T^2)/2+w;
                 
                if  abs( d1(k)-BBC)<ef/2  && alpha < vw*T  %%速度波门
                    track(j).assoi_point = [track(j).assoi_point;cc k];
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

%% 用第四次扫描的点迹继续判别航迹，采用椭圆波门 %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %第四次扫描的量测数
    tracknum = size(track,2);    %求得此时的航迹数
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %如果航迹在上一帧没有关联点迹，则对航迹进行直线外推两次，到d点
          X =  F*[track(j).seq(end,1) (track(j).seq(end,1)-track(j).seq(end-1,1))/T track(j).seq(end,2) (track(j).seq(end,2)-track(j).seq(end-1,2))/T]'; 
          P = F*[sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2]*F';
          X1 = F*X;
          P1 = F*P*F';
          Z1 = H*X1;
          S = H*P1*H'+R;
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
              if d(k) < 250 && alpha < thetamax
                  track(j).assoi_point = [track(j).assoi_point;d(k) k];
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
           P = F*[sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2]*F';
           X1 = F*X;
           P1 = F*P*F';%%预测协方差方程
           Z1 = H*X1;%%预测量测位置矩阵
           S = H*P1*H' +R;%%新息协方差矩阵
           for k = 1:samplenum
               data2 = cycle(i).sample(k,:)';
               d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
               a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
               b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
               c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
               alpha = acos((a + b - c)/(2*sqrt(a*b)));
               alpha = 180 - alpha * 180 / pi;
               if d(k) < 250 && alpha < thetamax %%9.21为参数查询X^2的分布表便可得到
                   track(j).assoi_point = [track(j).assoi_point;d(k) k];
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

%% 修正的hough变换 %%
Monte_Carlo=10;%Monte_Carlo仿真次数
target=5;%目标数
targets=5;%序列数
k=90;%sig分的个数
m=700;%p分的个数
T=1000;%雷达的量测距离

Np=1:k;
dNp=pi/k;%参数空间角度间隔
angle=(Np-1/2)*dNp;

dMp=2*T/m;%参数空间垂距间隔

success=zeros(Monte_Carlo,target);%目标航迹成功起始矩阵
fake(1:Monte_Carlo)=0;%目标航迹虚假起始矩阵
track_number(1:Monte_Carlo)=0;%总航迹起始数
A=zeros(k,m);%积累矩阵

 %计算相应的pn估值
 num=1;
    for i=1: tracknewnum
        L=size(tracknew(i).seq,1);
        for l=1:L
          for j=1:k
            P(num,j)=tracknew(i).seq(l,1)*cos(angle(j))+tracknew(i).seq(l,2)*sin(angle(j));
%             jiansuo(num)=i;
          end
          num=num+1;
        end
    end
    
  %对积累矩阵投票,看pn是否落在方格里，如果是就加一
    for i=1:k
        for j=1:m
            a=-T+(j-1)*dMp;
            b=-T+j*dMp;
%             a=(j-1)*dMp;
%             b=j*dMp;
           for h=1:size(P,1)
               if (P(h,i)>=a && P(h,i)<b) 
                   A(i,j)=A(i,j)+1;
               end
           end
        end
    end
    %% 模糊矩阵删选A矩阵 %%
% %%利用隶属度函数计算隶属度累积量
% [pm,~]=max(P);
% sita_m=90;%%最大误差范围
% 
% for i=1:k
%         for j=1:m
%             if  A(i,j)==0
%                 r=(j-1/2)*dMp;
%                 alph=(i-1/2)*dNp;
%             else
%                  alph=(i-1/2)*dNp;
%                  r=P(i,j);
%             end
%          zita_sita=;%%计算θ 方向上的均方差
%          zita_p=;%%计算 ρ方向上的均方差
%         if sita_m-(i-1)*dNp<3*zita_sita && pm-(j-1)*dMp<3*zita_p
%             A(i,j)= gaussmf(sita_m,[ zita_sita (i-1)*dNp])* gaussmf(pm,[zita_p (j-1)*dMp]);
%         else
%             A(i,j)=0;
%         end       
%         end
%  end
     %寻找投票数最大的五个参数
    for s=1:targets
        max=A(1,1);
        maxi=1;maxj=1;
        for i=1:k
            for j=1:m
                if A(i,j)>=max
                    max=A(i,j);
                    maxi=i;
                    maxj=j;
                end
            end
        end
        for i=maxi-5:maxi+5
            for j=maxj-5:maxj+5
                if i<=0
                    i=1;
                end
                if j<=0
                    j=1;
                end
                A(i,j)=0;
            end
        end
        P0(s)=-T+(maxj-1/2)*dMp;
%         x_init(s)=jiansuo(maxi);
        A0(s)=angle(maxi);
    end

%% 画图 %%
figure(1);
s = ['*', 's', '+', '.'];
hold on;
for i = 1:N
    plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
end

plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
     [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
for ii = 1:size(tracknew,2)
    data = tracknew(ii).seq;
    plot(data(:,1),data(:,2), '-');
end 
xlim([0, Xscope]);
ylim([0, Yscope]);
box on;

% figure(2);
% hold on;
% for ii = 1:size(tracknew,2)
%     data = tracknew(ii).seq;
%     plot(data(:,1),data(:,2), '-');
% end
% xlim([0, Xscope]);
% ylim([0, Yscope]);
% box on;
% 
% figure(3);
% hold on;
% plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
%      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
% for ii = 1:size(tracknew,2)
%     data = tracknew(ii).seq;
%     scatter(data(:,1),data(:,2),5,'r');
% end 
% for s=1:targets     
%         X=0:1:1000;
%         YS=(P0(s)-X*cos(A0(s)))/(sin(A0(s)));
%         plot(X,YS,'k');
%         hold on
% end
% xlim([0, Xscope]);
% ylim([0, Yscope]);
% box on;

% figure(3);
% s = ['*', 's', '+', '.'];
% hold on;
% for i = 1:N
%     plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
% end
% 
% plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
%      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
% for ii = 1:size(track1,2)
%     data = track1(ii).seq;
%     plot(data(:,1),data(:,2), '-');
% end 
% xlim([0, Xscope]);
% ylim([0, Yscope]);
% box on;
% 
% figure(4);
% s = ['*', 's', '+', '.'];
% hold on;
% for i = 1:N
%     plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
% end
% 
% plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
%      [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
% for ii = 1:size(track2,2)
%     data = track2(ii).seq;
%     plot(data(:,1),data(:,2), '-');
% end 
% xlim([0, Xscope]);
% ylim([0, Yscope]);
% box on;

% %% 航迹检测概率 %%
% syms x z
% f = exp(-x^2/2);
% 
% figure(3)
% hold on;
% sigma = 100:100:1000;
% m = length(sigma);
% for i =1:m
%     L = (vmin*T/(sqrt(2)*sigma(i)))^2;
%     H = (vmax*T/(sqrt(2)*sigma(i)))^2;
%     u = v*T/(sqrt(2)*sigma(i));
%     
%     f1 = 2/sqrt(pi)*int(f,0,sqrt((H-z^2)/2))*(exp(-(z-u)^2/2) + exp(-(z+u)^2/2));
%     f2 = 2/sqrt(pi)*int(f,0,sqrt((L-z^2)/2))*(exp(-(z-u)^2/2) + exp(-(z+u)^2/2));
%     f3 = (f)^(2*u^2);
%     p1 = 1/pi*(int(f1,0,sqrt(H)) - int(f2,0,sqrt(L)));
%     p2 = 2*u/sqrt(pi)*int(f3,0,tan(pi/4));
%     p(i) = p1*p2;
%     p(i) = double(p(i));
%     plot(sigma(i),p(i),'*');
% end
% % axis([100 1000 0 1]);
% xlim([100, 1000]);
% ylim([0, 10]);
% box on;
% i = 1:m;
% plot(sigma(i),p(i),'b*');
% xlable(100,1000);
% ylable(0,1);
% hold on


 
