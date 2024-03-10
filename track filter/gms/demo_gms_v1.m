% This is a demonstration of the adaptive GLMB filter with linear Gaussian 
% models. 


disp('Running adaptive GLMB filter with linear Gaussian models (v1) ...')
model = gen_model;
truth = gen_truth(model);
meas = gen_meas(model,truth);
est = run_filter_dpGLMB_v1(model,meas);
handles = plot_results(model,truth,meas,est);