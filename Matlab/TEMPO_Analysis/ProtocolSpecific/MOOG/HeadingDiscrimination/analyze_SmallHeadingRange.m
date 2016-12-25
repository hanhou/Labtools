aa = dlmread('HeadingDiscri_cum_vespsyneumetric.dat','',1,1);
dim = size(aa);

unique_heading1 = [-9,-3.5,-1.3,-0.5, 0, 0.5, 1.3, 3.5, 9]; % for monkey 1
unique_heading2 = [-16,-6.4,-2.6,-1, 0, 1, 2.6, 6.4, 16]; % for monkey 2

for i = 1 : dim(1)
    i
    if i <= 126
        unique_heading_big = unique_heading1;
        unique_heading_small = unique_heading1(2:end-1);
    else
        unique_heading_big = unique_heading2;
        unique_heading_small = unique_heading2(2:end-1);
    end
    xi_big = min(unique_heading_big) : 0.1 : max(unique_heading_big);
    xi_small = min(unique_heading_small) : 0.1 : max(unique_heading_small);
    
    % behavior 
	fit_data_psycho_big(1:9, 1) = unique_heading_big;  % full heading range
	fit_data_psycho_big(1:9, 2) = aa(i, 2 : 10);
	fit_data_psycho_big(1:9, 3) = aa(i,1); 
    
    fit_data_psycho_small(1:7, 1) = unique_heading_small;  % exclude the two biggest headings
	fit_data_psycho_small(1:7, 2) = aa(i, 2+1 : 10-1);
	fit_data_psycho_small(1:7, 3) = aa(i,1); 
    
    % neuron
    fit_data_neu_big(1:8, 1) = unique_heading_big(unique_heading_big~=0);  % full heading range
	fit_data_neu_big(1:8, 2) = aa(i, 2+9 : 10+9-1);
	fit_data_neu_big(1:8, 3) = aa(i,1); 
    
    fit_data_neu_small(1:6, 1) = unique_heading_small(unique_heading_small~=0);  % exclude the two biggest headings
	fit_data_neu_small(1:6, 2) = aa(i, 2+1+9 : 10-1+9-1);
	fit_data_neu_small(1:6, 3) = aa(i,1); 
    
    % fit curve now
    wichman_psy = pfit(fit_data_psycho_big,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
    Thresh_psy_big(i) = wichman_psy.params.est(2);
    psy_perf_big = [wichman_psy.params.est(1),wichman_psy.params.est(2)];

    wichman_psy = pfit(fit_data_psycho_small,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
    Thresh_psy_small(i) = wichman_psy.params.est(2);
    psy_perf_small = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
    
    wichman_neu = pfit(fit_data_neu_big,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
    Thresh_neu_big(i) = wichman_neu.params.est(2);
    neu_perf_big = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
    
    wichman_neu = pfit(fit_data_neu_small,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
    Thresh_neu_small(i) = wichman_neu.params.est(2); 
    neu_perf_small = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
    
    % calculate error, notice that only compare the middle points
    yi_psy_big = cum_gaussfit(psy_perf_big, xi_big);   
    %index_psy_big = [1,56,78,86,91,96,104,126,181];  
    index_psy_big = [56,78,86,91,96,104,126]; 
    err_psy_big(i) = sum( (yi_psy_big(index_psy_big) - fit_data_psycho_small(:,2)').^2 );
    
    yi_psy_small = cum_gaussfit(psy_perf_small, xi_small);
    index_psy_small = [1,23,31,36,41,49,71];
    err_psy_small(i) = sum( (yi_psy_big(index_psy_small) - fit_data_psycho_small(:,2)').^2 );
    
    yi_neu_big = cum_gaussfit(neu_perf_big, xi_big);
    heaing_big = unique_heading_big(unique_heading_big~=0);
    % index_neu_big = [1,56,78,86,96,104,126,181];
    index_neu_big = [56,78,86,96,104,126];
    err_neu_big(i) = sum( (yi_neu_big(index_neu_big) - fit_data_neu_small(:,2)').^2 );
    
    yi_neu_small = cum_gaussfit(neu_perf_small, xi_small);
    heaing_small = unique_heading_small(unique_heading_small~=0);
    index_neu_small = [1,23,31,41,49,71];
    err_neu_small(i) = sum( (yi_neu_small(index_neu_small) - fit_data_neu_small(:,2)').^2 );
       
end

% negative and positive infinite value means flat tuning
Thresh_neu_big(find(Thresh_neu_big<0|Thresh_neu_big>300))=300;
Thresh_neu_small(find(Thresh_neu_small<0|Thresh_neu_small>300))=300;

figure(1);
set(1,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
orient landscape;

subplot(2,2,1);
plot(Thresh_psy_small, Thresh_psy_big, 'o');
xlabel('narrow');
ylabel('full');
% xlim([0.1,5]);
% ylim([0.1,5]);
title('behavior threshold');

subplot(2,2,2);
loglog(Thresh_neu_small, Thresh_neu_big, 'o');
xlabel('narrow');
ylabel('full');
% xlim([0.3,300]);
% ylim([0.3,300]);
title('neuron threshold');

subplot(2,2,3);
plot(err_psy_small, err_psy_big, 'o');
xlabel('narrow');
ylabel('full');
% xlim([0.1,5]);
% ylim([0.1,5]);
title('behavior err');

subplot(2,2,4);
plot(err_neu_small, err_neu_big, 'o');
xlabel('narrow');
ylabel('full');
% xlim([0.3,300]);
% ylim([0.3,300]);
title('neuron err');