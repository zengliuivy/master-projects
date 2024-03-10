function data = simutrack(x0, y0, v, theta, sigma_ax, sigma_ay, sigma_r, sigma_theta, T, N)
%simutrack ��������ٶ��Ŷ�������ֱ���˶�Ŀ��Ķ�ά����.
%
%     'x0'          Ŀ����x�����ϵĳ�ʼλ��
%     'y0'          Ŀ����y�����ϵĳ�ʼλ��
%     'v'           Ŀ�������˶����ٶ�
%     'theta'       Ŀ���ٶ���x�᷽��ļнǣ���λ��
%     'sigma_ax'    x�᷽���������ٶ�
%     'sigma_ay'    y�᷽���������ٶ�
%     'sigma_r'     �������¾���Ĳ�����׼��
%     'sigma_theta' �������·�λ�Ĳ�����׼���λ��
%     'T'           �״�ɨ������
%     'N'           ��������
%
%     'data'        ����õ���N��Ŀ�꺽��
% ת��Ϊ����
theta = theta*pi/180;
sigma_theta = sigma_theta*pi/180;
% x��y�����ϵĳ�ʼ�ٶ�
vx0 = v * cos(theta);
vy0 = v * sin(theta);
% �Ŷ�Э�������
Q = [sigma_ax^2 0;0 sigma_ay^2];
Gamma = [T^2/2 0;T 0;0 T^2/2;0 T];
% ״̬ת�ƾ���
Phi = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
% ��������
H = [1 0 0 0;0 0 1 0];
% ������ʵ����
X(:,1) = [x0 vx0 y0 vy0]';
for m = 2:N
    X(:,m) = Phi*X(:,m-1)+Gamma*[sigma_ax*randn(1) sigma_ay*randn(1)]';
end
Pos = [X(1,:); X(3,:)];
% �������µ���ֵ
r0 = sqrt(X(1,:).^2+X(3,:).^2);
theta0 = atan(X(1,:)./X(3,:));  
% �Ӹ�˹����  
r = r0 + sigma_r*randn(1,N);  
theta = theta0 + sigma_theta*randn(1,N);    
% ������������������ת����ֱ������
x = r.*sin(theta)*exp(sigma_theta^2/2);
y = r.*cos(theta)*exp(sigma_theta^2/2);
data = [x', y'];
end

