function model= gen_model_line(T,P_D,P_S,lambda_c,clutter_P_S,clutter_P_D ,r_birth_clt,sigma_vel,sigma_e)

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

% % dynamical model parameters (CV model)
 model.T= T;                                           %sampling period采样周期
 model.A0= [ 1 model.T ; 0 1];                         %transition matrix一维转移矩阵                      
 model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];%二维转移矩阵
 model.B0= [ (model.T^2)/2; model.T ];%一维过程噪声分布矩阵
 model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];%二维过程噪声分布矩阵
 model.sigma_v = sigma_vel;%过程噪声标准差
 model.Q= (model.sigma_v)^2* model.B*model.B';         %process noise covariance过程噪声协方差
 model.pdvarfac_tg= 0.01^2;

% survival/death parameters
model.P_S=P_S;
model.Q_S= 1-model.P_S;

%这里不一样，这里设置了出生参数
% birth parameters
model.B_birth= diag([ 1000; 1000; 1000; 1000 ]);     % std of Gaussians
model.P_birth = model.B_birth*model.B_birth';      % cov of Gaussians
model.init_st = [90 , 10] ; %beta parameters of birth tracks


% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ sigma_e; sigma_e ]); 
model.D_r= diag([ 1*(pi/180); 5 ]); 
model.R= model.D*model.D';              %observation noise covariance

% detection parameters
model.P_D=P_D;             %to generate scenario, not known to the filter
model.Q_D= 1-model.P_D;     %probability of missed detection in measurements

%杂波参数设置不一样
% clutter parameters
model.range_c = [ -100000 100000; -100000 100000];                  %uniform clutter region
model.pdf_c = 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
model.lambda_cb = 1;                                         %birth rate for clutter targets
model.w_cb = 1;
model.clutter_P_S = clutter_P_S;                                     %survival probability for clutter targets
model.clutter_P_D = clutter_P_D ;                                    %detection probability for clutter targets
model.Lc_birth = 1; model.u_cb = 1; model.v_cb = 1;          %parameters of birth clutter
model.clutter_Nt = 300;                                      %number of clutter generators
model.lambda_c =lambda_c ;                                        %clutter rate to generate scenario, not known to the filter

%新的参数设置
% detection parameters
% model.P_D= .98;   %probability of detection in measurements
% model.Q_D= 1-model.P_D; %probability of missed detection in measurements
% model.clutter_P_D= .5;
% model.clutter_N_T= 20;
%model.lambda_c= model.clutter_N_T.*model.clutter_P_D;
model.range_r= [ -pi/2 pi/2; 0 200000 ];          %uniform clutter on r/theta
model.pdf_r= 1/prod(model.range_r(:,2)-model.range_r(:,1)); %uniform clutter density
model.B_birth_cb= diag([ 10000; 10000; 10000; 10000; 10*(pi/180) ]);                  %std of Gaussians
model.P_birth_cb= model.B_birth_cb*model.B_birth_cb';

model.P_S_clt= .90;
model.P_D_clt= .50;

model.T_clt= 1;
model.Q_clt= diag([1000 0 500 0 ].^2);%更改
model.pdvarfac_clt= 0.07^2;
% clutter parameters
model.D_clt = diag([20*(pi/180); 400]); %std for angle and range noise
model.R_clt = model.D_clt*model.D_clt';


model.T_birth_clt= 20;
model.L_birth_clt= zeros(model.T_birth_clt,1);
model.r_birth_clt= zeros(model.T_birth_clt,1);                                          %prob of birth for each LMB birth term
model.lambda_b_clt= cell(model.T_birth_clt,1);
model.m_birth_clt= cell(model.T_birth_clt,1);                                           %means of GM for each LMB birth term
model.B_birth_clt= cell(model.T_birth_clt,1);                                           %std of GM for each LMB birth term
model.P_birth_clt= cell(model.T_birth_clt,1);                                           %cov of GM for each LMB birth term
model.u_b_clt= zeros(model.T_birth_clt,1); model.v_b_clt= zeros(model.T_birth_clt,1);

for i=1:model.T_birth_clt
    model.L_birth_clt(i)= 1;
    model.r_birth_clt(i)= r_birth_clt;
    model.lambda_b_clt{i}(1)=1;
    model.m_birth_clt{i}(:,1)= [0; 0; 0; 0];%更改
    model.B_birth_clt{i}(:,:,1)= diag([200000; 0; 200000; 0]);%更改
    model.u_b_clt(i)= 5;
    model.v_b_clt(i)= 5;
    
    for g=1:model.L_birth_clt(i)
        model.P_birth_clt{i}(:,:,g)= model.B_birth_clt{i}(:,:,g)*model.B_birth_clt{i}(:,:,g)';
    end
end
