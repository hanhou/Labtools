%fit psychophysical data with gaussian

aa = dlmread('YongFitGau.txt');

fit_data_psycho(:, 1) = aa(:,1);  
fit_data_psycho(:, 2) = aa(:,2);
fit_data_psycho(:, 3) = 40; 

wichman_psy = pfit(fit_data_psycho,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
Thresh = wichman_psy.params.est(2);
Bias = wichman_psy.params.est(1);
psy_perf = [wichman_psy.params.est(1),wichman_psy.params.est(2)];

xi = min(aa(:,1)) : 0.1 : max(aa(:,1)); 
yi = cum_gaussfit(psy_perf, xi);

aa_out(:,1) = xi';
aa_out(:,2) = yi';
dlmwrite('YongFitGauM2_out.txt', aa_out); 
    
Thresh
Bias