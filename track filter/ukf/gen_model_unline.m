function model= gen_model_unline(T,P_D,P_S,lambda_c,clutter_P_S,clutter_P_D ,r_birth_clt,sigma_vel,sigma_turn)
% basic parameters
model.x_dim= 5;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= T;                         %sampling period
model.sigma_vel= sigma_vel;
model.sigma_turn= sigma_turn;   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim);
model.Q= model.B*model.B';
model.pdvarfac_tg= 0.01^2;
 
% survival/death parameters
model.P_S= P_S;
model.Q_S= 1-model.P_S;

% birth parameters (LMB birth model, single component only)
model.B_birth = diag([ 100; 100; 100; 100; 3*(pi/180) ]); 
model.P_birth =model. B_birth * model.B_birth' ; 

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 1*(pi/180); 5 ]);       %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% detection parameters
model.P_D= P_D;             %probability of detection in measurements
model.Q_D= 1-model.P_D;     %probability of missed detection in measurements
model.init_st = [90 , 10] ; %beta parameters of birth tracks

% clutter parameters
% model.lambda_c= 15;                                       %poisson average rate of uniform clutter (per scan)
model.range_c= [ -pi/2 pi/2; 0 20000 ];                      %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

model.lambda_cb = 1;                                %birth rate for clutter targets
model.w_cb = 1;
model.clutter_P_S = clutter_P_S;                            %survival probability for clutter targets
model.clutter_P_D = clutter_P_D ;                           %detection probability for clutter targets
model.Lc_birth = 1; model.u_cb = 1; model.v_cb = 1; %parameters of birth clutter
model.clutter_Nt = 300;                             %number of clutter generators
model.lambda_c= lambda_c;                                %clutter rate to generate scenario, not known to the filter

model.P_S_clt= .90;
model.P_D_clt= .50;

model.T_clt= 1;
model.Q_clt= diag([1000 0 500 0 0].^2);%
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
    model.m_birth_clt{i}(:,1)= [0; 0; 0; 0;0];
    model.B_birth_clt{i}(:,:,1)= diag([20000; 0; 20000; 0;0]);
    model.u_b_clt(i)= 5;
    model.v_b_clt(i)= 5;
    
    for g=1:model.L_birth_clt(i)
        model.P_birth_clt{i}(:,:,g)= model.B_birth_clt{i}(:,:,g)*model.B_birth_clt{i}(:,:,g)';
    end
end




