aa{1}=dlmread('metric1.txt');
aa{2}=dlmread('metric2.txt');
aa{3}=dlmread('metric3.txt');
aa{4}=dlmread('metric4.txt');
%aa{5}=dlmread('metric.txt');

for k =1 :4
    fit_data_psycho_cum{k}(:, 1) = aa{k}(:,1);  
    fit_data_psycho_cum{k}(:, 2) = aa{k}(:,2);  
    fit_data_psycho_cum{k}(:, 3) = aa{k}(:,3);     
      
    hold on;
    %calculate threshold
    wichman_psy = pfit(fit_data_psycho_cum{k},'plot_opt','plot without stats','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false'); 
    %wichman_psy = pfit(fit_data_psycho_cum,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false'); 
    threshold(k) =  wichman_psy.params.est(2) % This lists out the threshold values on the command window. ab feb-08.
end
% xmin=min([fit_data_psycho_cum{1}(:, 1),fit_data_psycho_cum{2}(:, 1),fit_data_psycho_cum{3}(:, 1),fit_data_psycho_cum{4}(:, 1)]);
% xmax=mmax([fit_data_psycho_cum{1}(:, 1),fit_data_psycho_cum{2}(:, 1),fit_data_psycho_cum{3}(:, 1),fit_data_psycho_cum{4}(:, 1)]);
xlim([-20 20]); %change this around to make different x axis limits. ab-jan-08
ylim([-0.2 1.2]); %extended y-axis to make more presentable. ab jan-08
axis square; %   set axes to be same length, while maintaining sphere shape of icons. ab jan-08

%     Thresh_psy(k) = wichman_psy.params.est(2);
%     Bias_psy(k) = wichman_psy.params.est(1);
%     psy_perf{k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
% 
% 
%     % fit curve
%     xi{k} = aa{k}(1,1) : 0.1 : aa{k}(end,1);   
%     beta = [0, 1.0];
%     yi{k} = cum_gaussfit(psy_perf{k}, xi{k});
% 
%     % plot figure   
%     plot(fit_data_psycho_cum{k}(:, 1), fit_data_psycho_cum(:, 2), 'bo');
%     hold on;
%     plot(xi, yi, 'r-');
%     ylim([0 1]);

