%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningPlot_GenData.m -- Generates model or noisy model data
%      sets based on models used for Curvefit.
%-- pwatkins, 5/04
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_GenData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, ...
    StartOffset, StopOffset, PATH, FILE);

% pwatkins - START data generation options

Curvefit_defines;

unique_stim_type = [CURVEFIT_OCCULAR_STIM];  % which stim types to generate data for

% fitting_models = [@Curvefit_head_centered_ng @Curvefit_eye_centered_ng];
% fitting_bounds = [@Curvefit_head_centered_bounds_ng @Curvefit_eye_centered_bounds_ng];
% fitting_tuning_models = [@Curvefit_cos_tuning_7p @Curvefit_cos_tuning_7p];
fitting_models = [ @Curvefit_eye_centered_5gf];
fitting_bounds = [ @Curvefit_eye_centered_bounds_5gf];
fitting_tuning_models = [ @Curvefit_cos_tuning_5p];
% fitting_models = [@Curvefit_head_centered_ng ];
% fitting_bounds = [@Curvefit_head_centered_bounds_ng ];
% fitting_tuning_models = [@Curvefit_cos_tuning_7p ];

% parameters to use with fitting model to generate data.
% each row corresponds to each fitting model in fitting_models array.
% empty matrix or row with all -1's means to use x0 guess from fitting_bounds.
fitting_model_params = [0 0 3 100/pi 25 20/pi 20];  

num_unique_elevation = 23;  % number of latitudes to generate (including poles)
num_unique_azimuth = 44;    % number of longitutdes to generate
num_data_sets = 1;   % how many data sets to create
num_complete_repititions = 5;  % how many repititions per data set (per point)
unique_gaze_angle = [-20 0 20];  % which gaze angles to generate data for
num_gazes = length(unique_gaze_angle);
fixation_types = [CURVEFIT_VARY_FIXATION_X];  % which fixation types to generate data for.
                                              % set to empty array for random per data set.

% START - for fitting see DirectionTuningPlot_Curvefit backdoor options
    %backdoor_dir = fullfile('work','tempo_backdoor','gendata');
    backdoor_dir = fullfile('work','tempo_backdoor');
    % create the backdoor directory, incase it does not exist.
    % for some stupid reason matlab can not create an absolute directory.
    tmp = pwd;
    cd(matlabroot);
    [s,mess,messid] = mkdir(backdoor_dir);
    cd(tmp);
    timestamp = datestr(clock,30);
    backdoor_dir = fullfile(matlabroot, backdoor_dir);
    backdoor_save_pre = ['gendata_' timestamp '_'];
    backdoor_load_ext = '_load.mat';
    backdoor_batch_file = fullfile(backdoor_dir, ['GenData_' timestamp '.m']);
% END - for fitting see DirectionTuningPlot_Curvefit backdoor options

% END data generation options

% START generate data

% create the data points about the sphere
unique_azimuth = (0:1:num_unique_azimuth-1)'/num_unique_azimuth*360;
unique_elevation = (-num_unique_elevation+1:2:num_unique_elevation-1)'/(num_unique_elevation-1)*90;
[el az] = meshgrid(unique_elevation(2:end-1), unique_azimuth);  % exclude the poles
unique_point_azimuth = [0 az(:)' 0];  % put the poles back
unique_point_elevation = [-90 el(:)' 90];
num_unique_points = length(unique_point_elevation);  % or azimuth

unique_point_azimuth_r = unique_point_azimuth/180*pi;
unique_point_elevation_r = unique_point_elevation/180*pi;
unique_gaze_angle_r = unique_gaze_angle/180*pi;
unique_gaze_angle_r = -unique_gaze_angle_r;   % use same convention as MOOG experiments

% pass this into fitting function.
% we are generating data, so we do not care about SSE.
fake_stim_data = zeros(num_gazes, num_unique_points, num_complete_repititions+1);

% allocate space for generated data.
% see DirectionTuningPlot_Curvefit for an explanation of the data structure.
resp_mat = zeros(max(unique_stim_type), num_gazes, num_unique_points, num_complete_repititions+1);

% initialize global variable kludges that have to do with fitting.
curvefit_gaze_fits = zeros(num_gazes,num_unique_points);
curvefit_gaze_residuals = zeros(num_gazes,num_unique_points);
curvefit_gaze_sse = zeros(1,num_gazes);
curvefit_gaze_mean_sse = zeros(1,num_gazes);

% create a dummy batch file.
fh_batch = fopen(backdoor_batch_file, 'w');

% create specified number of data sets for each fitting type.
% each data set could contain data for multiple stim types.
num_stims = length(unique_stim_type);
num_fits = length(fitting_models);
for k=1:num_data_sets
    for j=1:num_fits
        
        % randomize gaze direction, if it was not specified.
        % only allow one gaze direction for all stim types.
        if length(fixation_types) == 0
            if rand(1,1) > 0.5
                fixation_type = CURVEFIT_VARY_FIXATION_X;
            else
                fixation_type = CURVEFIT_VARY_FIXATION_Y;
            end
        else
            fixation_type = fixation_types(1);
        end
        
        for i=1:num_stims

            % if no parameters were supplied, use the random ones
            % generated by the bounding function.  Do not randomize
            % for each repitition, since these data are expected to
            % have variance accounted for by noise only.
            R = feval(fitting_bounds(j), [], num_gazes);
            if length(fitting_model_params) == 0 | ...
                    fitting_model_params(j,:) == -ones(1,size(fitting_model_params,2))
                current_model_params = R{CURVEFIT_BOUNDS_X0};
            else
                current_model_params = fitting_model_params(j,1:length(R{CURVEFIT_BOUNDS_X0}));
            end
            
            % generate data over all repitions
            for m=1:num_complete_repititions
                feval( fitting_models(j), current_model_params, fitting_tuning_models(j), ...
                    unique_point_azimuth_r, unique_point_elevation_r, fake_stim_data, ...
                    num_complete_repititions, unique_gaze_angle_r, fixation_type );
                
                resp_mat(unique_stim_type(i), :, :, m+1) = add_noise(curvefit_gaze_fits);
            end
            
            % store the mean in the first repitition index, as usual.
            resp_mat(unique_stim_type(i), :, :, 1) = ...
                mean(resp_mat(unique_stim_type(i), :, :, 2:num_complete_repititions+1),4);
        end
        
        % store the generated data for this data set.
        % variables stored need to agree with those in DirectionTuningPlot_Curvefit
        FILE = sprintf('%s_%s_%s_%05d', backdoor_save_pre, func2str(fitting_models(j)), ...
            gendata_noise_type, k);
        backdoor_load_file = fullfile(backdoor_dir, [FILE backdoor_load_ext]);
        save(backdoor_load_file, 'resp_mat', 'unique_stim_type', 'unique_gaze_angle_r', ...
            'unique_point_azimuth_r', 'unique_point_elevation_r', 'num_complete_repititions', ...
            'num_unique_points', 'fixation_type', 'num_gazes', 'unique_azimuth', ...
            'unique_elevation', 'unique_point_azimuth', 'unique_point_elevation', ...
            'unique_gaze_angle');
        sprintf('data saved to %s', backdoor_load_file)
        
        % write the data set to the batch file
        fprintf( fh_batch, 'Z:\\nada %s ''Output Data for Curve Fitting'' -1 -1 4 500 5 -500 1 1\n', FILE );
    end
end

fclose(fh_batch);
return;


%-----------------------------------------------------------------------------------------------------------------------
%-- add_noise -- Internal function to add noise to ideal data.
%-----------------------------------------------------------------------------------------------------------------------
function noisy_data = add_noise( trial_data )

Curvefit_defines;

gendata_noise_type = 'no_noise';
noisy_data = trial_data;

return;


%-----------------------------------------------------------------------------------------------------------------------
% END OF LINE
%-----------------------------------------------------------------------------------------------------------------------
