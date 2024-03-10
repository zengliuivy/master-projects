%����Hough�任
close all
clear all

fig_num =0;
% ɨ�������ɨ������
N = 3;
T = 5; %��

% �����ǵ������η�������
Xscope = 10^5;
Yscope = 10^5;

% Ŀ���˶�����
v = 500;     % 500m/s
theta =0;   % ˮƽ��x���˶�

sigmax=1; %x�᷽���������ٶ�
sigmay=0.6; %y�᷽���������ٶ�

%%%%%��������%%%%%%%%%%%
S=v/5;%%������������
w=sqrt(S);%%���þ�����

% ����۲��׼���뷽λ�ǹ۲��׼��
sigma_r = w;%%��˹������׼��
sigma_theta = 0.3;

Vmin = 2*v/3;
Vmax = 3*v/2;

% �����ǵ������η��������ڵ��Ӳ�ƽ����
renbuda =120;

%���ؿ������ʵ�����
mc=10;
for cy=1:mc
% k=90;%sig�ֵĸ���
% m=500;%p�ֵĸ���
% ָ������ɨ����Ӳ�������ÿ�����ڵ���Ŀ���Ӳ��ɷֲ����ֲ��ľ�ֵ�������С
% �Լ���λ������Ӳ����ĳ˻�ȷ��
K = poissrnd(renbuda, 1, N);

% cycle1=[];
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
% radar5 = simutrack(55000, 15000, v, theta, 0, 0, 0, 0, T, N);

 radar1 = simutrack(55000, 55000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N); %4��2��
radar2 = simutrack(45000, 45000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar3 = simutrack(35000, 35000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar4 = simutrack(45000, 25000, v, theta, sigmax, sigmay, sigma_r, sigma_theta, T, N);
radar5 = simutrack(55000, 15000, v, theta,sigmax, sigmay, sigma_r, sigma_theta, T, N);
% tic
% %�����
%% ÿ��ɨ�����õ㼣����sample�е�ǰ5���㱻�趨ΪĿ��� %%
i = 0;
for k = K
    i = i + 1;
     cycle(i).sample = [rand(k,1)*Xscope rand(k,1)*Yscope];       %cycleΪ�ṹ��   �洢�Ӳ�
    cycle(i).sample = [radar1(i,:); radar2(i,:); radar3(i,:);
        radar4(i,:); radar5(i,:); cycle(i).sample];
end
% for l=1:N
%  cycle1=[cycle1;cycle(l).sample];
% end
Z1= cycle(1).sample;
Z2= cycle(2).sample;
Z3= cycle(3).sample;
% Z4= cycle(4).sample;

%����ʽ����
number=1;
for i=1:size(Z1,1)
    for j=1:size(Z2,1)
        for h=1:size(Z3,1)
            x11=Z1(i,1);y11=Z1(i,2);
            x21=Z2(j,1);y21=Z2(j,2);
            x31=Z3(h,1);y31=Z3(h,2);
            V12=sqrt((x11-x21)*(x11-x21)+(y11-y21)*(y11-y21))/T;
            V23=sqrt((x31-x21)*(x31-x21)+(y31-y21)*(y31-y21))/T;
            if Vmin<=V12 && V12<=Vmax && Vmin<=V23 && V23<=Vmax
                Zc1(number,1:2)=[x11,y11];            
                Zc2(number,1:2)=[x21,y21];
                Zc3(number,1:2)=[x31,y31];
                number=number+1;
            end
        end
    end
end

k=200;%%sig�ֵĸ���
m=400;%p�ֵĸ���
Np=1:k;
dNp=pi/k;%�����ռ�Ƕȼ��
% dNp=180/k;%�����ռ�Ƕȼ��
angle=(Np-1/2)*dNp;

dMp=2*T/m;%�����ռ䴹����
%���㴹��
for i=1:number-1
    x11=Zc1(i,1);y11=Zc1(i,2);
    x21=Zc2(i,1);y21=Zc2(i,2);
    x31=Zc3(i,1);y31=Zc3(i,2);
    for j=1:k 
        P1(i,j)=x11*cos(angle(j))+y11*sin(angle(j));
        P2(i,j)=x21*cos(angle(j))+y21*sin(angle(j));
        P3(i,j)=x31*cos(angle(j))+y31*sin(angle(j));
    end
end
count=1;
%�����ֺ���,Ѱ�ҹ����
for i=1:number-1
    for j=1:k
        deta_P1(i,j)=P2(i,j)- P1(i,j);
        deta_P2(i,j)=P3(i,j)- P2(i,j);
        if j>1
        if (deta_P1(i,j)/deta_P1(i,j-1))<0 
            deta_01(count,:)=[i,j];
             count=count+1;
        end
        end
    end
end

%Ѱ�ҷ���Ҫ��������
num=1;
for i=1:count-1
    for j=2:k
         deta_02=deta_P2(deta_01(i,1),j);
         deta_12=deta_P2(deta_01(i,1),j-1);
         if (deta_02/ deta_12)<0 
             %���������
                x11=Zc1(deta_01(i,1),1);y11=Zc1(deta_01(i,1),2);
                x21=Zc2(deta_01(i,1),1);y21=Zc2(deta_01(i,1),2);
                x31=Zc3(deta_01(i,1),1);y31=Zc3(deta_01(i,1),2);
                slope12=deta_P1(deta_01(i,1),deta_01(i,2));
                slope23=deta_02;
             if abs(j-deta_01(i,2))<=30 && (slope23/slope12)>0 %%%�о�
                  para(num,:)=[x11,y11,x21,y21,x31,y31];
                  num=num+1;
             end
         end
    end
end
%��ͼ
%% hough�任���
fig_num = fig_num+1;
figure(fig_num);
% s = ['*', 's', '+'];
% hold on;
% for i = 1:N
%     plot(cycle(i).sample(6:end,1), cycle(i).sample(6:end,2), s(i));
% end
plot([radar1(:,1); radar2(:,1); radar3(:,1); radar4(:,1); radar5(:,1)],...
     [radar1(:,2); radar2(:,2); radar3(:,2); radar4(:,2); radar5(:,2)], 'o');
 hold on;

for h=1:num-1
    X=[para(h,1),para(h,3),para(h,5)];
    YS=[para(h,2),para(h,4),para(h,6)];
    plot(X,YS,'b');
    hold on
end
% 
%   toc
% disp(['����ʱ��: ',num2str(toc)]);
end

