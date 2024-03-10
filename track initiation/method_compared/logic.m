% %% Track Start  ���Ŀ�꺽����ʼ  3/4�߼�

clear all; close all;
warning('off')
fig_num =0;
% ɨ�������ɨ������
N = 4;
M = 3;
T = 5; %��

% �����ǵ������η�������
Xscope = 10^5;
Yscope = 10^5;

% Ŀ���˶�����
v = 500;     % 500m/s
theta = 0;   % ˮƽ��x���˶�

sigmax=1; %x�᷽���������ٶ�
sigmay=0.6; %y�᷽���������ٶ�

%%%%%��������%%%%%%%%%%%
S=v/5;%%������������
 w=sqrt(S);%%���þ�����

% ����۲��׼���뷽λ�ǹ۲��׼��
sigma_r = w;%%��˹������׼��
sigma_theta = 0.3;

% �����ǵ������η��������ڵ��Ӳ�ƽ����
renbuda =120;
    
% ָ������ɨ����Ӳ�������ÿ�����ڵ���Ŀ���Ӳ��ɷֲ����ֲ��ľ�ֵ�������С
% �Լ���λ������Ӳ����ĳ˻�ȷ��
K = poissrnd(renbuda, 1, N);

% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% ���ƹ��������е��������С�ٶȡ������ٶȺ���С���ٶȣ���������ɨ��ļн�
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
thetamax =30;
yama=25;
% thetamax = pi/2;

%���ⷽ��
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
% R = [sigma_r^2  sigma_theta^2;sigma_theta^2 sigma_r^2];%%%���������������й�
R = [25000  0;0 25000];
P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
%   P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];

%���ؿ������ʵ�����

mc=1;

for h=1:mc
    tic
%�����
Temp(1).sample=[]; Temp(2).sample=[];Temp(3).sample=[];Temp(4).sample=[];%%%�ṹ�壬���ڴ洢�Ӳ����˺������


%% �������5��Ŀ��ĺ���(��������) %%
% radar1 = simutrack(55000, 55000, v, theta, 0, 0, sigma_r, sigma_theta, T, N); %4��2��
% radar2 = simutrack(45000, 45000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar3 = simutrack(35000, 35000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar4 = simutrack(45000, 25000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar5 = simutrack(55000, 15000, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar1 = simutrack(55000, 55000, v, theta, 0, 0, 0, 0, T, N); %4��2��
% radar2 = simutrack(45000, 45000, v, theta, 0, 0, 0, 0, T, N);
% radar3 = simutrack(35000, 35000, v, theta, 0, 0, 0, 0, T, N);
% radar4 = simutrack(45000, 25000, v, theta, 0, 0, 0, 0, T, N);
% radar5= simutrack(55000, 15000, v, theta, 0, 0, 0, 0, T, N);
% 
 radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4��2��
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
%% ÿ��ɨ�����õ㼣����sample�е�ǰ5���㱻�趨ΪĿ��� %%
i = 0;
for k = K
    i = i + 1;
    cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycleΪ�ṹ��   �洢�Ӳ���
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
end

%% �õ�һ��ɨ��ĵ㼣������ʱ���� %%
for i = 1:size(cycle(1).sample, 1)
    track(i).seq = cycle(1).sample(i,:);
%     track(i).shouldadd = [];
    track(k).assoi_point = [];      %�洢�뺽�������ĵ㼣�Ĺ���ֵ
end

%% �õڶ���ɨ��ĵ㽨�����ܺ��� %%
for i = 2
    tracknum = size(track,2);      %�����̬������
    tracknum_temp = tracknum;
    samplenum = size(cycle(i).sample,1);     %��õڶ�֡������㼣��
    
    D = zeros(tracknum,samplenum);      %�洢��̬����������Ĺ���ֵ
    %% ���㱾��ɨ������е㼣����̬�����Ĺ���ֵ %%
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
        
        %% ����̬����δ�����ĵ㼣��Ϊ�µ���̬����ͷ
        if flag == 0
            tracknum_temp =tracknum_temp + 1;
            
            track(tracknum_temp).seq = cycle(i).sample(j,:);
            track(tracknum_temp).assoi_point = [];
        end
    end
    
    %% �ɹ����㼣�б𣬶���̬�������д��� %%
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
    
    %% ���ຽ�� %%
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

%% �õ���֡ɨ��ĵ㼣�����ɿ����� %%
for i = 3
    tracknum = size(track,2);
    samplenum = size(cycle(i).sample,1);
    for j = 1:tracknum
        track(j).assoi_point = [];
        num = size(track(j).seq,1);     %����ÿ�������еĵ㼣��
        %% ��̬������ֻ��һ���㼣ʱ %% 
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
        %% ��̬�����ж���һ���㼣ʱ %%
        if num > 1
            data = track(j).seq(end-1,:);     %�����еĵ����ڶ����㼣
            data1 = track(j).seq(end,:);      %�����е����һ���㼣
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
     
    %% �Ժ����������� %% 
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
    
%     %% ���ڿɿ������еĵ㼣�е��������ģ���Ϊ�ɹ���ʼ���º��� %%
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

%% �õ��Ĵ�ɨ��ĵ㼣�����б𺽼� %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %���Ĵ�ɨ���������
    tracknum = size(track,2);    %��ô�ʱ�ĺ�����
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %�����������һ֡û�й����㼣����Ժ�������ֱ������
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
    
%     %% ���ຽ�� %%
%     track3 = [];
%     track3num = 0;
%     for j = 1:tracknum
%         if ~isempty(track(j).seq)
%             track3num = track3num + 1;
%             
%             track3(track3num).seq = track(j).seq;
%         end
%     end
    %%  ���ڿɿ������еĵ㼣�е��������ģ���Ϊ�ɹ���ʼ���º��� %%   
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
% disp(['����ʱ��: ',num2str(toc)]);

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
