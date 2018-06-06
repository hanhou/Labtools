% population analysis for 3D models
% from data loading -> analysis
% LBY 20170607


models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
FnameCode = 2;
Protocol = 1;


popul_model_load_data(models,FnameCode, Protocol); % load data
popul_modelFit_evalu_analysis(models,Protocol) % evaluate the models