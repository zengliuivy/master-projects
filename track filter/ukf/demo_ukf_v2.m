% This is a demonstration of the adaptive GLMB filter with non-linear Gaussian 
% models.

disp('Running adaptive GLMB filter with non-linear Gaussian models (v2) ...')
model = gen_model;
truth = gen_truth(model);
meas = gen_meas(model,truth);
est = run_filter_v2(model,meas);
handles = plot_results(model,truth,meas,est);