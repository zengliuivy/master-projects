clear all; close all;
warning('off')

fig_num =0;
%���ؿ���������
mc=1;

 for h=1:mc

%        tic
% %�����
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
renbuda =90;

% ָ������ɨ����Ӳ�������ÿ�����ڵ���Ŀ���Ӳ��ɷֲ����ֲ��ľ�ֵ�������С
% �Լ���λ������Ӳ����ĳ˻�ȷ��
K = poissrnd(renbuda, 1, N);

% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% ���ƹ��������е��������С�ٶȡ������ٶȺ���С���ٶȣ���������ɨ��ļн�
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
vw=6;%%ת����ٶȣ���λΪ��/s��ת����ٶȵ������벨������й�
thetamax =30;
yama=25;
% thetamax = pi/2;

%���ⷽ��
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
 %P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
 R = [10000  0;0 10000];%%%���������������й�
% R = [sigma_r^4  sigma_theta^4;sigma_theta^4 sigma_r^4];
P0=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];


Temp(1).sample=[]; Temp(2).sample=[];Temp(3).sample=[];Temp(4).sample=[];%%%�ṹ�壬���ڴ洢�Ӳ����˺������
if ~exist('track','var')
else
for r=1:size(track,2)
track(r).assoi_point=[];track(r).seq=[];
end
end

%% �������5��Ŀ��ĺ���(��������) %%

radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4��2��
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
%% ÿ��ɨ�����õ㼣����sample�е�ǰ5���㱻�趨ΪĿ��� %%
i = 0;
for k = K
    i = i + 1;
     cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycleΪ�ṹ��   �洢�Ӳ�
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
    
end

%% �õ�һ��ɨ��ĵ㼣������ʱ���� %%
for i = 1:size(cycle(1).sample, 1)
    track(i).seq = cycle(1).sample(i,:); %%%�洢���еĵ�һ��ɨ��Ľ��
%     track(i).shouldadd = [];
    track(k).assoi_point = [];      %�洢�뺽�������ĵ㼣�Ĺ���ֵ
end

%% �õڶ���ɨ��ĵ㽨�����ܺ��� %%
for i = 2
    tracknum = size(track,2);      %�����̬������������һ��ɨ�����������������
    tracknum_temp = tracknum;
    samplenum = size(cycle(i).sample,1);     %��õڶ���ɨ�����е�����㼣��
    
    D = zeros(tracknum,samplenum);      %�洢��̬����������Ĺ���ֵ
    %% ���㱾��ɨ������е㼣����̬�����Ĺ���ֵ %%
    for j = 1:samplenum
        data = cycle(i).sample(j,:); %%��ǰ�ڶ���ɨ�������ֵ
        for k = 1:tracknum
            if size(track(k).seq,1) > 0 
                data1 = track(k).seq; %%��ǰ��һ��ɨ�������ֵ
                D(k,j) = (data(1)-data1(1))^2 + (data(2)-data1(2))^2; %%����֮��ľ���ƽ��
                alpha(k,j)=AngX([data(1)-data1(1);data(2)-data1(2)]);%%%Ŀ�귽�����
            end                 
        end
    end
    
    for j = 1:samplenum
        flag = 0;
        
        for k = 1:tracknum
            if D(k,j) >= (vmin*T-w)^2 && D(k,j) <= (vmax*T+w)^2 && alpha(k,j)<90%%�������,����Ŀ���˶��������ƣ�����x������
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
%         for k = 1:tracknum_temp
        L = size(track(k).assoi_point,1);
        if L == 1 %%����Ѿ�����
            j = track(k).assoi_point(end,2);%�ڶ���ɨ��ĵڼ����㼣
            track(k).seq = [track(k).seq;cycle(i).sample(j,:)];%%�������ϵĵ㼣���ݶ�����һ��seq��
        end
        if L > 1%����кü����������ѡ��������Ǹ���
            min = track(k).assoi_point(1,:);
            for j = 2:L
                if (track(k).assoi_point(j,1) - (v*T)^2) < (min(1) - (v*T)^2)
                    min = track(k).assoi_point(j,:);
                end
            end
            track(k).seq = [track(k).seq;cycle(i).sample(min(2),:)];
        end
        if L == 0 %���û������ֱ�����
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

%% �õ���֡ɨ��ĵ㼣�����������������β��� %%
for i = 3
    tracknum = size(track,2);
    samplenum = size(cycle(i).sample,1);
    for j = 1:tracknum
        track(j).assoi_point = [];
        num = size(track(j).seq,1);     %����ÿ�������еĵ㼣��
        %% ��̬������ֻ��һ���㼣ʱ���ڵ����������������ظ���������������ʹ��Բ�β���%% 
        if num == 1   
            data = track(j).seq(end,:);
            for k = 1:samplenum
                data1 = cycle(i).sample(k,:);
                d(k) = (data1(1) - data(1)*v*T )^2 + (data1(2) - data(2))^2;%%%%��Ŀ�������������
                alpha(k,j)=AngX([data1(1)-data(1)*v*T ;data1(2)-data(2)]);%%%Ŀ�귽�����
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
        %% ��̬�����ж���һ���㼣ʱ %%
        if num > 1
            data = track(j).seq(end-1,:);     %�����еĵ����ڶ����㼣
            data1 = track(j).seq(end,:);      %�����е����һ���㼣
            %%����
%             X1 = [data1(1) (data1(1)-data(1))/T data1(2) (data1(2)-data(2))/T]';%%�ú��������� x���ٶ� ������ Y���ٶ�
% 
%             Z1 = H*F*X1;%%%���Ƶ�����
            %�������˲��Ľ�
            X=[data(1);v;data(2);0];
%         P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
   P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
            x_ =F*X;          %Ԥ�⵱ǰ״̬
            P1 = F*P*F';
            % ���²�
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data1(1);data1(2)] - H*x_);
            P= (1-k*H)*P1;
           
            x_ = F*x;            %Ԥ�⵱ǰ״̬
            Z1=H*x_;
            P1 = F*P*F';
            S = H*P1*H'+ R;%%��ϢЭ�������
            for k = 1:samplenum 
                data2 = cycle(i).sample(k,:)';%%%%����Ȧ��������
%                 d1(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
%                 a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
%                 b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
%                 c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
%                 alpha = acos((a + b - c)/(2*sqrt(a*b)));
%                 alpha = 180 - alpha * 180 / pi;
                d1(k) =sqrt((data2(1)-data1(1))^2 + (data2(2)-data1(2))^2);%�������⵽b�����
                BBC=sqrt((Z1(1)-data1(1))^2 + (Z1(2)-data1(2))^2);%����bc�����
                bc=[data2(1)-data1(1);data2(2) - data1(2)];%%bc'����
                bbc=[Z1(1) - data1(1);Z1(2) - data1(2)];%%bc����
               cc=sqrt((Z1(1)-data2(1))^2+(Z1(2)-data2(2))^2);%cc'֮��ľ���
%                 %�������ת���
                alpha=abs(AngX(bc)-AngX(bbc));
%                a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
%                b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
%                c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
%                alpha = acos((a + b - c)/(2*sqrt(a*b)));
%                alpha = 180 - alpha * 180 / pi;
%            plot([data(1);data1(1);Z1(1)],[data(2);data1(2);Z1(2)],'o');hold on;
%            plot(radar1(:,1),radar1(:,2));
                %������������뾶
                ef=(vmax-vmin)*T+((amax-amin)*T^2)/2+w;
                 
                if  abs( d1(k)-BBC)<ef/2 && alpha < vw*T  %%�ٶȲ���  
                    track(j).assoi_point = [track(j).assoi_point;cc k];
                      Temp(3).sample=[Temp(3).sample;cycle(i).sample(k,:)];
                end
            end
            L = size(track(j).assoi_point,1);
            if L == 1
                k = track(j).assoi_point(end,2);
                track(j).seq = [track(j).seq;cycle(i).sample(k,:)];%%����µ�����㵽����
            end
            if L > 1 %%����ֹһ�����������
                min = track(j).assoi_point(1,:);
                for k = 2:L
                    if track(j).assoi_point(k,1) < min(1)
                        min = track(j).assoi_point(k,:);
                    end
                end
                track(j).seq = [track(j).seq;cycle(i).sample(min(2),:)];%%%ѡ��������Ƶ�c��������������
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
    
end

%% �õ��Ĵ�ɨ��ĵ㼣�����б𺽼���������Բ���� %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %���Ĵ�ɨ���������
    tracknum = size(track,2);    %��ô�ʱ�ĺ�����
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %�����������һ֡û�й����㼣����Ժ���������������
          X =  F*[track(j).seq(end,1) (track(j).seq(end,1)-track(j).seq(end-1,1))/T track(j).seq(end,2) (track(j).seq(end,2)-track(j).seq(end-1,2))/T]'; 
          data1 = H*X;
%           track(j).seq = [track(j).seq;data1'];
           X=[track(j).seq(end,1);v;track(j).seq(end,2);0];
%            P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
%  R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
%     R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
%     R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
           P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            x_ =F*X;          %Ԥ�����֡״̬
            P1 = F*P*F';
            % ���²�
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data1(1);data1(2)] - H*x_);
            P= (1-k*H)*P1;
                  
             x_ = F*x; %Ԥ�����֡���Ƶ�
            Z1=H*x_;
        P1 = F*P*F';
           S = H*P1*H' +R;%%��ϢЭ�������
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
%            P1 = F*P*F';%%Ԥ��Э�����
%            Z1 = H*X1;%%Ԥ������λ�þ���
           
              %�������˲��Ľ�
            X=[track(j).seq(end-2,1);v;track(j).seq(end-2,2);0];
           P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
     %       P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            x_ =F*X;          %Ԥ��ڶ�֡״̬
            P1 = F*P*F';
            % ���²�
            k = P1*H'/(H*P1*H'+R);
            x = x_ + k*([data(1);data(2)] - H*x_);
            P= (1-k*H)*P1;
           
            x_ = F*x;            %Ԥ�����֡״̬
            P1 = F*P*F';
            %����
             k = P1*H'/(H*P1*H'+R);
             x = x_ + k*([data1(1);data1(2)] - H*x_);
             P= (1-k*H)*P1;
                  
             x_ = F*x; %Ԥ�����֡���Ƶ�
            Z1=H*x_;
            P1 = F*P*F';
           S = H*P1*H' +R;%%��ϢЭ�������
%     X=[data(1);v;data(2);0];
%          P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
% % P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%��ʼЭ�������
% %     R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
% %     R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
% %     R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];
%             x_ =F*X;          %Ԥ�����֡״̬
%             P1 = F*P*F';
%             % ���²�
%             k = P1*H'/(H*P1*H'+R);
%             x = x_ + k*([data1(1);data1(2)] - H*x_);
%             P= (1-k*H)*P1;
%            
%             x_ = F*x;            %Ԥ�⵱ǰ״̬
%             Z1=H*x_;
%             S = H*P1*H'+ R;%%��ϢЭ�������
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
               if d(k) < yama && alpha < thetamax %%������ѯX^2�ķֲ����ɵõ�
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

    tracknew = [];
    tracknewnum = 0;
    for j = 1:tracknum
        if size(track(j).seq,1) > 2
            tracknewnum = tracknewnum + 1;
            tracknew(tracknewnum).seq = track(j).seq;
        end
    end
end

%% ģ���ĵ�hough�任 %%
target=5;%Ŀ����
k=200;%sig�ֵĸ���
m=400;%p�ֵĸ���
Ln=150000;%�״���������

Np=1:k;
dNp=pi/k;%�����ռ�Ƕȼ��
% dNp=180/k;%�����ռ�Ƕȼ��
angle=(Np-1/2)*dNp;

dMp=2*Ln/m;%�����ռ䴹����
[~,Seq]=size(tracknew);
for i=1:Seq
A(i).seq=zeros(k,m);%���۾����ʼ������Ϊȫ�����
A1(i).seq=zeros(k,m);%���۾����ʼ������Ϊȫ�����
end
AA=zeros(k,m);
 
%% hough�任 %%
  %������Ӧ��pn��ֵ
    for i=1: Seq
        L=size(tracknew(i).seq,1);
        for l=1:L
          for j=1:k
            Pn(i).sample(l,j)=tracknew(i).seq(l,1)*cos(angle(j))+tracknew(i).seq(l,2)*sin(angle(j));
          end
        end
    end

    %�Ի��۾���ͶƱ,��pn�Ƿ����ڷ��������Ǿͼ�һ��������ֵɸѡ����ҪΪ�˼��ټ�����
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
   
   %�Ի��۾���ͶƱ,����ģ�������Ⱥ���
%   sita_m=deg2rad(90);%%�����Χ
  sita_m=1;%%�����Χ
  pm=1400;%%%p�������Χ
%   zita_sita=1;
%    zita_p=200;
  index=1;
 for l=1: Seq
      for i=1:k
        for j=1:m
            if A1(l).seq(i,j)>2%%%���ýϵ͵���ֵ��Ϊ�˼��ټ�����
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
                 if abs(Pn(l).sample(h,g)-a)<=2*dMp  %%�Ƿ����뷽�񸽽��ĸ�����
                         pii=angle(g);
                         pj=Pn(l).sample(h,g);%���÷�ֵ����Ϊ�Ա���  
                         for ii=1:(2*sita_m/dNp+1)
                            sita_ijl(ii)=pii-sita_m+(ii-1)*dNp;
                         end
                         for jj=1:(2*pm/dMp+1)
                             p_ijl(jj)=pj-pm+(jj-1)*dMp;
                         end
                          zita_sita=std(sita_ijl,0);%%����� �����ϵľ�����
                          zita_p=std(p_ijl,0);%%���� �ѷ����ϵľ�����
                           uk=uk+gaussmf(pii,[ zita_sita b])* gaussmf(pj,[zita_p a]);%%%�����Ⱥ���
%                         A(l).coordinate(i,j)=index;%%��������
%                         Co(index).seq(num,:)=[tracknew(l).seq(h,1),h,i,j];%%%�����洢x���꣬���ڻ�ԭ�켣
%                         num=num+1;
                 end
                   end
               end
                  A(l).seq(i,j)=A(l).seq(i,j)+uk;%%%%�ۻ�����ļ���
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
%     xlabel('��');
%     ylabel('��');
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
 %% ��ֵ��ȷ�Ϻ���%%
 
 num=1;track_hough=[];
 %%Ѱ��ȫ����ֵ
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
             if max_sum> max_mean*0.68%%������ֵ
%            if max_sum>8%%������ֵ
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


%% ��ͼ %%
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


% % % %% hough�任���
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
% disp(['����ʱ��: ',num2str(toc)]);
 end
% % 



