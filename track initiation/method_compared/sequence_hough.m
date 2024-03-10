%����Hough�任
clear all; close all;
warning('off')

fig_num =0;
% ɨ�������ɨ������
N = 12;
M = 3;
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

% �����ǵ������η��������ڵ��Ӳ�ƽ����
renbuda =30;

%���ؿ������ʵ�����
mc=1;

for h=1:mc
% ָ������ɨ����Ӳ�������ÿ�����ڵ���Ŀ���Ӳ��ɷֲ����ֲ��ľ�ֵ�������С
% �Լ���λ������Ӳ����ĳ˻�ȷ��
K = poissrnd(renbuda, 1, N);

cycle1=[];
%% �������5��Ŀ��ĺ���(��������) %%

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
for l=1:N
 cycle1=[cycle1;cycle(l).sample];
end

target=5;%Ŀ����
targets=6;%������
k=200;%sig�ֵĸ���
m=400;%p�ֵĸ���
L=150000;%�״��������
Pd=1;%������

Np=1:k;
dNp=pi/k;%�����ռ�Ƕȼ��
angle=(Np-1/2)*dNp;

dMp=2*L/m;%�����ռ䴹����

A=zeros(k,m);%���۾���
 coordinate=zeros(k,m);
% for monte=1:Monte_Carlo
ii=1;max_sum=10;
% while max_sum>targets
for s=1:targets
  %��������� Hough�任
  P=[];
    for i=1:size(cycle1,1)
        for j=1:k
            P(i,j)=cycle1(i,1)*cos(angle(j))+cycle1(i,2)*sin(angle(j));
        end
    end

 
    %�Ի��۾���ͶƱ
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
                   coordinate(i,j)=index;%%��������
                   Co(index).seq(num,:)=[cycle1(h,1),h,i,j];%%%�����洢x���꣬���ڻ�ԭ�켣
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
     cycle1(all(cycle1==0,2),:)=[];%%%ɾ������ֵ����Ӧ����
 end
    
    %��ͼ
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
%     fprintf('�ܹ������ĺ�������Ϊ %f\n',s); 
%   toc
% disp(['����ʱ��: ',num2str(toc)]);  
end

        