clear all; close all;
warning('off')

% ɨ�������ɨ������
N = 4;
M = 3;
T = 5; %��

% �����ǵ������η�������
Xscope = 10^3;
Yscope = 10^3;

% Ŀ���˶�����
v = 5;     % 5m/s
theta = 0;   % ˮƽ��x���˶�

sigmax=0.5; %x�᷽���������ٶ�
sigmay=0.06; %y�᷽���������ٶ�

%%%%%��������%%%%%%%%%%%
S=v/5;%%������������
w=sqrt(S);%%���þ�����

% ����۲��׼���뷽λ�ǹ۲��׼��
sigma_r =w;%%��˹������׼��
sigma_theta = 0.3;

% �����ǵ������η��������ڵ��Ӳ�ƽ����
renbuda = 30;
    
% ָ������ɨ����Ӳ�������ÿ�����ڵ���Ŀ���Ӳ��ɷֲ����ֲ��ľ�ֵ�������С
% �Լ���λ������Ӳ����ĳ˻�ȷ��
K = poissrnd(renbuda, 1, N);


% normrnd(0,w)*Xscope  normrnd(0,w)*Yscope
% ���ƹ��������е��������С�ٶȡ������ٶȺ���С���ٶȣ���������ɨ��ļн�
vmin = 2*v/3;
vmax = 3*v/2;
amax = (vmax)/T;
amin=0;
vw=6;%%ת����ٶȣ���λΪ��/s
thetamax = 30;
% thetamax = pi/2;

%���ⷽ��
H = [1 0 0 0;0 0 1 0];
F = [1 T 0 0; 0 1 0 0;  0 0 1 T;0 0 0 1];
R = [160 0;0 250];


%% �������5��Ŀ��ĺ���(��������) %%
radar1 = simutrack(550, 550, v, theta, 0, 0, sigma_r, sigma_theta, T, N); %4��2��
radar2 = simutrack(450, 450, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar3 = simutrack(350, 350, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar4 = simutrack(450, 250, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
radar5 = simutrack(550, 150, v, theta, 0, 0, sigma_r, sigma_theta, T, N);
% radar1 = simutrack(550, 550, v, theta, 0, 0, 0, 0, T, N); %4��2��
% radar2 = simutrack(450, 450, v, theta, 0, 0, 0, 0, T, N);
% radar3 = simutrack(350, 350, v, theta, 0, 0, 0, 0, T, N);
% radar4 = simutrack(450, 250, v, theta, 0, 0, 0, 0, T, N);
% radar5 = simutrack(550, 150, v, theta, 0, 0, 0, 0, T, N);

%  radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4��2��
% radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
% radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
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
                d(k) = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
                alpha(k,j)=AngX([data1(1)-data(1);data1(2)-data(2)]);%%%Ŀ�귽�����
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
        %% ��̬�����ж���һ���㼣ʱ %%
        if num > 1
            data = track(j).seq(end-1,:);     %�����еĵ����ڶ����㼣
            data1 = track(j).seq(end,:);      %�����е����һ���㼣
            %���Ƽ���
            X = [data1(1) (data1(1)-data(1))/T data1(2) (data1(2)-data(2))/T]';%%�ú��������� x���ٶ� ������ Y���ٶ�
            P = [sigmax^2 sigmax^2/T 0 0;sigmax^2/T 2*sigmax^2/T^2 0 0;0 0 sigmay^2 sigmay^2/T;0 0 sigmay^2/T sigmay^2/T^2];
            X1 = F*X;
            P1 = F*P*F';
            Z1 = H*X1;%%%���Ƶ�����


            S = H*P1*H' + R;%%��ϢЭ�������
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
                %�������ת���
                alpha=abs(AngX(bc)-AngX(bbc))*2;
            
                %������������뾶
                ef=(vmax-vmin)*T+((amax-amin)*T^2)/2+w;
                 
                if  abs( d1(k)-BBC)<ef/2  && alpha < vw*T  %%�ٶȲ���
                    track(j).assoi_point = [track(j).assoi_point;cc k];
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

%% �õ��Ĵ�ɨ��ĵ㼣�����б𺽼���������Բ���� %%
for i = 4
    samplenum = size(cycle(i).sample,1);       %���Ĵ�ɨ���������
    tracknum = size(track,2);    %��ô�ʱ�ĺ�����
    for j = 1:tracknum
        L = size(track(j).assoi_point,1);
        data = track(j).seq(end,:);
        if L == 0    %�����������һ֡û�й����㼣����Ժ�������ֱ���������Σ���d��
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
           P1 = F*P*F';%%Ԥ��Э�����
           Z1 = H*X1;%%Ԥ������λ�þ���
           S = H*P1*H' +R;%%��ϢЭ�������
           for k = 1:samplenum
               data2 = cycle(i).sample(k,:)';
               d(k) = (data2 - Z1)'*inv(S)*(data2 - Z1);
               a = (data2(1) - data1(1))^2 + (data2(2) - data1(2))^2;
               b = (data1(1) - data(1))^2 + (data1(2) - data(2))^2;
               c = (data2(1) - data(1))^2 + (data2(2) - data(2))^2;
               alpha = acos((a + b - c)/(2*sqrt(a*b)));
               alpha = 180 - alpha * 180 / pi;
               if d(k) < 250 && alpha < thetamax %%9.21Ϊ������ѯX^2�ķֲ����ɵõ�
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

%% ������hough�任 %%
Monte_Carlo=10;%Monte_Carlo�������
target=5;%Ŀ����
targets=5;%������
k=90;%sig�ֵĸ���
m=700;%p�ֵĸ���
T=1000;%�״���������

Np=1:k;
dNp=pi/k;%�����ռ�Ƕȼ��
angle=(Np-1/2)*dNp;

dMp=2*T/m;%�����ռ䴹����

success=zeros(Monte_Carlo,target);%Ŀ�꺽���ɹ���ʼ����
fake(1:Monte_Carlo)=0;%Ŀ�꺽�������ʼ����
track_number(1:Monte_Carlo)=0;%�ܺ�����ʼ��
A=zeros(k,m);%���۾���

 %������Ӧ��pn��ֵ
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
    
  %�Ի��۾���ͶƱ,��pn�Ƿ����ڷ��������Ǿͼ�һ
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
    %% ģ������ɾѡA���� %%
% %%���������Ⱥ��������������ۻ���
% [pm,~]=max(P);
% sita_m=90;%%�����Χ
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
%          zita_sita=;%%����� �����ϵľ�����
%          zita_p=;%%���� �ѷ����ϵľ�����
%         if sita_m-(i-1)*dNp<3*zita_sita && pm-(j-1)*dMp<3*zita_p
%             A(i,j)= gaussmf(sita_m,[ zita_sita (i-1)*dNp])* gaussmf(pm,[zita_p (j-1)*dMp]);
%         else
%             A(i,j)=0;
%         end       
%         end
%  end
     %Ѱ��ͶƱ�������������
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

%% ��ͼ %%
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

% %% ���������� %%
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


 
