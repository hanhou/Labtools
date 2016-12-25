
%num_hidden_units_trials = {[50:10:150] [175:25:300] [350:50:500]};
%num_hidden_units_trials = {[50:150]};
%num_hidden_units_trials = {[150]};

%num_unique_azimuth_trials = [10:2:22];
%num_unique_elevation_trials = [6:12];
num_unique_azimuth_trials = [8];
num_unique_elevation_trials = [5];

%num_hidden_units_reps = [50];
%num_hidden_units_reps = [2];

total_cnt = 1;
ind = 1;

global num_hidden_units;  % this is a hack to set this param from a top level looping script
global num_hidden_units_load;
global num_unique_elevation;  % this is a hack to set this param from a top level looping script
global num_unique_elevation_load;
global num_unique_azimuth;  % this is a hack to set this param from a top level looping script
global num_unique_azimuth_load;

fig_rows = 4;
fig_cols = 3;

for mm=1:length(num_unique_azimuth_trials)
    for ii=1:length(num_hidden_units_trials)
        cur_trials = num_hidden_units_trials{ii};
        for jj=1:length(cur_trials)

            %num_pts(total_cnt) = ...
            %    (num_unique_azimuth_trials(mm))*(num_unique_elevation_trials(mm)-2)+2;
            num_pts(total_cnt) = cur_trials(jj);
            
            sel = ind:ind+num_hidden_units_reps(ii)-1;
            
            figure(775);
            subplot(fig_rows,fig_cols,total_cnt);
            data = ang_prefcv_gt_ang_prefc_cnt(sel)./cur_trials(jj);
            hist(data); 
            title(sprintf('%d hidden units - mean %.3f, stdv %.3f', ...
                cur_trials(jj), mean(data), std(data)));
            all_mean_ang(total_cnt) = mean(data);
            all_std_ang(total_cnt) = std(data);

            figure(776);
            subplot(fig_rows,fig_cols,total_cnt);
            data = network_mse_all(sel);
            hist(data); 
            title(sprintf('%d hidden units - mean %.2e, stdv %.2e', ...
                cur_trials(jj), mean(data), std(data)));
            all_mean_mse(total_cnt) = mean(data);
            all_std_mse(total_cnt) = std(data);
            
            figure(777);
            subplot(fig_rows,fig_cols,total_cnt);
            data = pref_visual_mag_gt_pref_vestib_mag_cnt(sel)./cur_trials(jj);
            hist(data); 
            title(sprintf('%d hidden units - mean %.3f, stdv %.3f', ...
                cur_trials(jj), mean(data), std(data)));
            all_mean_pref(total_cnt) = mean(data);
            all_std_pref(total_cnt) = std(data);
            
            total_cnt = total_cnt+1;
            ind = ind+num_hidden_units_reps(ii);
            
        end
    end
end
