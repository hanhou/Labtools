
%num_hidden_units_trials = {[50:10:150] [175:25:300] [350:50:500]};
%num_hidden_units_trials = {[50:150]};
num_hidden_units_trials = {[130]};

%num_unique_azimuth_trials = [10:2:22];
%num_unique_elevation_trials = [6:12];
num_unique_azimuth_trials = [8];
num_unique_elevation_trials = [5];

num_hidden_units_reps = [50];
%num_hidden_units_reps = [2];

total_cnt = 1;
ind = 1;
data_dir = 'data';

global num_hidden_units;  % this is a hack to set this param from a top level looping script
global num_hidden_units_load;
global num_unique_elevation;  % this is a hack to set this param from a top level looping script
global num_unique_elevation_load;
global num_unique_azimuth;  % this is a hack to set this param from a top level looping script
global num_unique_azimuth_load;

for mm=1:length(num_unique_azimuth_trials)
    for ii=1:length(num_hidden_units_trials)
        cur_trials = num_hidden_units_trials{ii};
        for jj=1:length(cur_trials)
            for kk=1:num_hidden_units_reps(ii)
                num_hidden_units = 0;
                num_hidden_units_load = cur_trials(jj);
                num_unique_azimuth = 0;
                num_unique_azimuth_load = num_unique_azimuth_trials(mm);
                num_unique_elevation = 0;
                num_unique_elevation_load = num_unique_elevation_trials(mm);
                run_network_model;
                
                pref_visual_mag_gt_pref_vestib_mag_cnt(total_cnt) = pref_visual_mag_gt_pref_vestib_mag;
                ang_prefcv_gt_ang_prefc_cnt(total_cnt) = ang_prefcv_gt_ang_prefc;
                network_mse_all(total_cnt) = network_mse;
                network_norm_mse_all(total_cnt) = network_norm_mse;
                network_msereg_all(total_cnt) = network_msereg;
                network_norm_msereg_all(total_cnt) = network_norm_msereg;
                ang_pref_all((total_cnt-1)*num_hidden_units+1:total_cnt*num_hidden_units) = ang_pref;
                total_cnt = total_cnt + 1;
                
                save_file = sprintf('save_hidden%d_az%d_el%d_trial%d.mat', num_hidden_units_load, ...
                    num_unique_azimuth_load, num_unique_elevation_load, kk);
                save(fullfile(data_dir, save_file), 'training_inputs', 'training_outputs', 'net', ...
                    'norm_rng', 'num_hidden_units', 'range_az', 'range_el', 'range_gaze_az', ...
                    'range_gaze_el', 'num_neurons', 'num_eye_neurons', 'P', 'P_eye', ...
                    'minp', 'maxp', 'mint', 'maxt', 'num_unique_azimuth', 'num_unique_elevation');
            end
        end
    end
end

save(fullfile(data_dir,'loop_save.mat'),'num_hidden_units_trials','num_hidden_units_reps', ...
    'num_unique_azimuth_trials', 'num_unique_elevation_trials', 'ang_pref_all', ...
    'pref_visual_mag_gt_pref_vestib_mag_cnt','ang_prefcv_gt_ang_prefc_cnt',...
    'network_mse_all','network_norm_mse_all','network_msereg_all','network_norm_msereg_all');
