function data = simutrack(x0, y0, v, theta, sigma_ax, sigma_ay, sigma_r, sigma_theta, T, N)
%simutrack 仿真带加速度扰动的匀速直线运动目标的二维航迹.
%
%     'x0'          目标在x方向上的初始位置
%     'y0'          目标在y方向上的初始位置
%     'v'           目标匀速运动的速度
%     'theta'       目标速度与x轴方向的夹角，单位度
%     'sigma_ax'    x轴方向的随机加速度
%     'sigma_ay'    y轴方向的随机加速度
%     'sigma_r'     极坐标下距离的测量标准差
%     'sigma_theta' 极坐标下方位的测量标准差，单位度
%     'T'           雷达扫描周期
%     'N'           采样点数
%
%     'data'        仿真得到的N点目标航迹
% 转变为弧度
theta = theta*pi/180;
sigma_theta = sigma_theta*pi/180;
% x和y方向上的初始速度
vx0 = v * cos(theta);
vy0 = v * sin(theta);
% 扰动协方差矩阵
Q = [sigma_ax^2 0;0 sigma_ay^2];
Gamma = [T^2/2 0;T 0;0 T^2/2;0 T];
% 状态转移矩阵
Phi = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
% 测量矩阵
H = [1 0 0 0;0 0 1 0];
% 构造真实航迹
X(:,1) = [x0 vx0 y0 vy0]';
for m = 2:N
    X(:,m) = Phi*X(:,m-1)+Gamma*[sigma_ax*randn(1) sigma_ay*randn(1)]';
end
Pos = [X(1,:); X(3,:)];
% 极坐标下的数值
r0 = sqrt(X(1,:).^2+X(3,:).^2);
theta0 = atan(X(1,:)./X(3,:));  
% 加高斯噪声  
r = r0 + sigma_r*randn(1,N);  
theta = theta0 + sigma_theta*randn(1,N);    
% 将加噪声的数据重新转换到直角坐标
x = r.*sin(theta)*exp(sigma_theta^2/2);
y = r.*cos(theta)*exp(sigma_theta^2/2);
data = [x', y'];
end

