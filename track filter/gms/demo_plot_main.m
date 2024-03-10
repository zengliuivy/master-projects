

% This code is based on the code framework written by Professor BA-NGU VO: 
% http://ba-ngu.vo-au.com/publications.html

% the _common folder needs to be added to the MATLAB path
% MC_max_runs determines the Monte Carlo runs

%=============================
% Scenario  

MC_max_runs =5;  % Monte Carlo runs
Outliers = 0;       % NO outliers, Outliers = 0, otherwise, Outliers = 1
lamda = 30;  % Poisson average rate of uniform clutter (per scan)
pd = 0.95;   % Detection probability

%=================================================
% Parameters of environment
ospa_c = 100;                                     % OSPA parameters
track_time=100;                                   % tracking time
filter_max_number =1;                            % the number of fitlers
T=1; 
win_len= 10;
cutoff = 100;

%=======================================================
 disp('Running adaptive GLMB filter with non-linear Gaussian models (v1) ...')
% start Monte Carlo runs
for MC_runs=1:MC_max_runs
                                                      
   
                                                      
    
    % generate model
  
   model = gen_model;
   truth = gen_truth(model);
   meas = gen_meas(model,truth);
    %==========================================================
    %                    start filtering
    %==========================================================
    
    % Prepare memory
    for huancun_i = 1:filter_max_number
        save_est(huancun_i). d1 = zeros(1,track_time); % OSPA errors
        save_est(huancun_i). d2 = cell(1,length(win_len)); % OSPA errors
        save_est(huancun_i). n = zeros(1,track_time); % the number of targets erros
        save_est(huancun_i). p = zeros(1,track_time); % ºÏ≤‚¬ 
        save_est(huancun_i). l = zeros(1,track_time); % ‘”≤®
    end
    
    % filtering the data
      filter_n=3;
      est =  run_filter_v2(model,meas); 
      [ospa,ospa2]=plot_results_MC(model,truth,meas,est);
  
     
        % save OSPA data
        save_est(filter_n).d1=ospa;
        save_est(filter_n).d2=ospa2;
        save_est(filter_n).n=est.N' ;
     %   save_est(filter_n).pd=est.pd;
     %   save_est(filter_n).cutter=est.lambda_c;
        
    % save data for plot section
     p_i =filter_n;
        save_data(p_i).d1=save_data(p_i).d1 + save_est(p_i).d1;
        save_data(p_i).n=save_data(p_i).n + save_est(p_i).n;
        for i = 1:length(win_len)
        save_data(p_i).d2{i}=save_data(p_i).d2{i} + save_est(p_i).d2{i};
        end
     %   save_data(p_i).pd=save_data(p_i).pd + save_est(p_i).pd;
      %  save_data(p_i).cutter=save_data(p_i).cutter + save_est(p_i).cutter;
    
end
%%save
%save save_data.mat save_data
