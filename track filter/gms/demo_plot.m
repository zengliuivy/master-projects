
% clear all
% clc,close all

% This code is based on the code framework written by Professor BA-NGU VO: 
% http://ba-ngu.vo-au.com/publications.html

% the _common folder needs to be added to the MATLAB path
% MC_max_runs determines the Monte Carlo runs

%=============================
% Scenario  

MC_max_runs = 1;  % Monte Carlo runs
Outliers = 0;       % NO outliers, Outliers = 0, otherwise, Outliers = 1
lamda = 30;  % Poisson average rate of uniform clutter (per scan)
pd = 0.95;   % Detection probability

%=================================================
% Parameters of environment
ospa_c = 100;                                     % OSPA parameters
track_time=100;                                   % tracking time
filter_max_number = 1;                            % the number of fitlers
T=1; 
win_len= 10;
cutoff = 100;

for i=1:filter_max_number
    plot_data(i).d1=zeros(track_time,3);           % OSPA errors
    plot_data(i).d2 =cell(1,length(win_len)); % OSPA errors
    plot_data(i).d2{1}=zeros(3,track_time); 
    plot_data(i).n=zeros(1,track_time);           % the number of targets erros
    plot_data(i).t=zeros(1,track_time);           % time cost
end

%=======================================================
 disp('Running adaptive GLMB filter with linear Gaussian models (v1) ...')
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
        save_est(huancun_i). p = zeros(1,track_time); % ¼ì²âÂÊ
        save_est(huancun_i). l = zeros(1,track_time); % ÔÓ²¨
    end
    
    % filtering the data
    for filter_n= 1:filter_max_number 
           
        % filtering
        if filter_n == 1
            est =  run_filter_v1(model,meas); % AMDP-J-GLMB
            est3 = est;
        elseif filter_n == 2
            est = run_filter_dpGLMB_v1(model,meas); % DP-GLMB
            est4 = est;
        elseif filter_n == 3
            est =run_filterGLMB(model,meas); % GLMB
        end
       [ospa,ospa2]=plot_results_MC(model,truth,meas,est);
  
     
        % save OSPA data
        save_est(filter_n).d1=ospa;
        save_est(filter_n).d2=ospa2;
        save_est(filter_n).n=est.N' ;
        
    end
   
    % save data for plot section
    for p_i = 1:filter_n
        plot_data(p_i).d1=plot_data(p_i).d1 + save_est(p_i).d1;
        plot_data(p_i).n=plot_data(p_i).n + save_est(p_i).n;
        for i = 1:length(win_len)
        plot_data(p_i).d2{i}=plot_data(p_i).d2{i} + save_est(p_i).d2{i};
        end
    end
end
%%save
save plot_data.mat plot_data

%==========================================================
%                    plot figures(MC_max_runs+1)
%==========================================================

for i=1:filter_max_number
    plot_data(i).d1 = plot_data(i).d1 /(MC_max_runs);
    plot_data(i).n = plot_data(i).n /( MC_max_runs);
     for j = 1:length(win_len)
     plot_data(i).d2{j} = plot_data(i).d2{j} /(MC_max_runs);
     end
end

%-----------------------------------------------------

figure(1)
c = {'-b.','-ro','-g>','-k+'};
number = truth.N;
plot(number),hold on
for i = 1:filter_max_number
    plot(plot_data(i).n,c{i})
end
legend('True number','AMDP-J-GLMB','DP-GLMB','Robudt-GLMB')
xlabel('Time (s)'),ylabel('Number of targets')
axis auto;

figure(2)
c={'b','r','g','k','m','c'};
for i = 2:filter_max_number
    subplot(3,1,1); plot(1:meas.K,plot_data(i).d1(:,1),c{i}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
    hold on
    subplot(3,1,2); plot(1:meas.K,plot_data(i).d1(:,2),c{i}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
    hold on
    subplot(3,1,3); plot(1:meas.K,plot_data(i).d1(:,3),c{i}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
    hold on
end
xlabel('Time (s)'),ylabel('OSPA, p=2, c=100')
legend('AMDP-J-GLMB','DP-GLMB','GLMB')

figure(3)
c={'-b.','-r.','-g','-k+','-m','-c'};
for j = 1:filter_max_number
    subplot(3,1,1);
    for i = 1:length(win_len)
      plot(1:truth.K,plot_data(j).d2{i}(1,:),c{j}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA$^{(2)}$ Dist','interpreter','latex');
       windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];
    end
    legend(windowlengthlabels,'interpreter','latex');hold on;

   subplot(3,1,2);
   for i = 1:length(win_len)
    plot(1:truth.K,plot_data(j).d2{i}(2,:),c{j}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA$^{(2)}$ Loc','interpreter','latex');
    windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];hold on;
   end

   subplot(3,1,3);
   for i = 1:length(win_len)
    plot(1:truth.K,plot_data(j).d2{i}(3,:),c{j}); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA$^{(2)}$ Card','interpreter','latex');
    windowlengthlabels{i} = ['$L_w = ' int2str(win_len(i)) '$'];hold on;
   end
end
xlabel('Time (s)','interpreter','latex')
legend('AMDP-J-GLMB','DP-GLMB','GLMB')

% average performance per scan per run
c = {'AMDP-J-GLMB','DP-GLMB','GLMB'};
for i=1:filter_max_number
    the_OSPA(i) = sum(plot_data(i).d1(:,1))/size(plot_data(i).d1,2);
    disp([c{i},' : ','     OSPA is ',num2str(the_OSPA(i),4)])
end
 a=plot_data(3).n ./truth.N';
  b=plot_data(2).n ./truth.N';
  
sum=0;
for i=1:100
    sum=sum+b(i);
end