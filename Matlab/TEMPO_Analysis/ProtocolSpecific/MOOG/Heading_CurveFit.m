%-----------------------------------------------------------------------------------------------------------------------
%-- Heading_CurveFit.m -- Fits response to model data as a function of azimuth and elevation for MOOG 3D tuning expt
%--	GCD, 6/27/03
%-- Modified for Vary_Fixation protocol   CRF + Yong, 12/03
%-- Modified for curve fitting - pwatkins, 4/04
%---Modified for fitting by Zack, 6/02/04
%-----------------------------------------------------------------------------------------------------------------------
function Heading_CurveFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, ...
    StartOffset, StopOffset, PATH, FILE, OPT_FILE);

Curvefit_defines;
script_name = 'DirectionTuningPlot_Curvefit';

% pwatkins - load an options file to specify curve fitting specific
% options.  this is an optional argument (not built into tempo_gui or
% batch_gui)
pause off;
if nargin < 13
    % use a default options file in the same directory as the script
    OPT_FILE = which(script_name);
    % get the directory that the script is in
    OPT_FILE = regexprep(OPT_FILE, ['\' filesep '[^\' filesep ']*$'], '');
    % use a default options file name
    OPT_FILE = fullfile( OPT_FILE, 'Curvefit_opt_def.m' );
    %OPT_FILE = fullfile( 'Z:', filesep, 'LabTools', filesep, 'Matlab', filesep, 'TEMPO_Analysis', filesep, 'ProtocolSpecific', filesep, 'MOOG', filesep, 'Curvefit_opt_def.m' );
else
    pause off;
end
% read the whole options file and eval into workspace
opt_file = textread(OPT_FILE, '%s', 'whitespace', '', 'bufsize', CURVEFIT_MAX_OPT_FILE_SIZE);
eval(cell2mat(opt_file), 'error([''OPT_FILE eval failed in '' script_name '' with error '' lasterr])');
%sprintf('successfully loaded options file %s', OPT_FILE)

% START load data
if ~run_load, return; end

if ~use_backdoor_load  % there is no backdoor, you are not reading this.

    Path_Defs;
    ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

    %get the column of values for azimuth and elevation and stim_type
    temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
    temp_elevation = data.moog_params(ELEVATION,:,MOOG);
    temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
    temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
    temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
    temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);

    %deal with data files where FIX_X and FIX_Y not defined
    temp_fix_x(isnan(temp_fix_x)) = 0;
    temp_fix_y(isnan(temp_fix_y)) = 0;

    %now, get the firing rates for all the trials
    temp_spike_rates = data.spike_rates(SpikeChan, :);

    %get indices of any NULL conditions (for measuring spontaneous activity
    null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

    %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
    trials = 1:length(temp_azimuth);		% a vector of trial indices
    bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
    if ( bad_trials ~= NaN)
        select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
    else
        select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
    end

    azimuth = temp_azimuth(~null_trials & select_trials);
    elevation = temp_elevation(~null_trials & select_trials);
    stim_type = temp_stim_type(~null_trials & select_trials);
    amplitude = temp_amplitude(~null_trials & select_trials);
    fix_x     = temp_fix_x(~null_trials & select_trials);
    fix_y     = temp_fix_y(~null_trials & select_trials);
    spike_rates = temp_spike_rates(~null_trials & select_trials);

    unique_azimuth  = munique(azimuth');
    unique_elevation = munique(elevation');
    unique_stim_type = munique(stim_type');
    unique_amplitude = munique(amplitude');
    unique_fix_x    =  munique(fix_x');
    unique_fix_y    =  munique(fix_y');

    % kludge to get the direction of the gaze angles
    if length(unique_fix_y) == 1
        gaze_angle = fix_x;
        temp_gaze_angle = temp_fix_x;
        fixation_type = CURVEFIT_VARY_FIXATION_X;
    else
        gaze_angle = fix_y;
        temp_gaze_angle = temp_fix_y;
        fixation_type = CURVEFIT_VARY_FIXATION_Y;
    end
    unique_gaze_angle = munique(gaze_angle');

    % kludge to get the unique points about the sphere
    num_unique_points = 0;
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) );
            if any(select)
                num_unique_points = num_unique_points + 1;
                unique_point_azimuth(num_unique_points) = unique_azimuth(i);
                unique_point_elevation(num_unique_points) = unique_elevation(j);
            end
        end
    end

    % create a working data structure that contains the mean and repitition
    % data.  data structure is indexed by (in order of indices):
    %   index type      valid indices          indexed by
    %   -------------------------------------------------
    %   stim type       [1,3]                  stim_type
    %   gaze angle      [1,num_gazes]          index into unique_gaze_angle
    %   point number    [1,num_unique_points]  index into unique_point_azimuth
    %                                                 and unique_point_elevation
    %   repitition      [1,num_reps+1]         repitition number
    %
    % Because there could be incomplete repititions over all the saved trials,
    % keep a counter of the minimum number of complete repitions.  Index 1 into
    % the repitition dimension stores the mean response over all repitions.
    % The purpose of the point arrays is so we only save the points around
    % the sphere, and do not create empty points for all azimuth angles at
    % the poles.
    resp_mat = [];
    num_complete_repititions = CURVEFIT_MAX_REPITITIONS;
    num_gazes = length(unique_gaze_angle);
    max_stim_type = max(unique_stim_type);
    for n=1:length(unique_stim_type)
        for k=1:num_gazes
            for i=1:num_unique_points
                select = logical( ...
                    (azimuth==unique_point_azimuth(i)) & (elevation==unique_point_elevation(i)) & ...
                    (gaze_angle==unique_gaze_angle(k)) & (stim_type==unique_stim_type(n)) );
                if any(select)
                    % store the mean data over all repititions.
                    % calculate the mean after we know number of complete
                    % repititions over all points.
                    %resp_mat(unique_stim_type(n), k, i, 1) = mean(spike_rates(select));

                    % remember the minimum number of complete repititions
                    num_reps = sum(select);
                    if num_complete_repititions > num_reps
                        num_complete_repititions = num_reps;
                    end

                    % store the trial data for each repitition
                    %DEBUGGING
%                    unique_stim_type(n)
%                    k
%                    i
%                    num_reps
%                    select
                    
                    
                    resp_mat(unique_stim_type(n), k, i, 2:num_reps+1) = spike_rates(select);
                else
                    % this really shouldn't happen for a good data set
                    num_complete_repititions = 0;
                end
            end
        end
    end

    % only keep complete repitions in resp_mat.
    % NOTE: for consistency, always index resp_mat or matrices derived from
    % resp_mat with the size of each dimension, instead of simply :,:,:
    resp_mat = resp_mat(1:max_stim_type, 1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1);

    % now calculate the mean, based only on the complete repititions.
    resp_mat(1:max_stim_type, 1:num_gazes, 1:num_unique_points, 1) = ...
        mean(resp_mat(1:max_stim_type, 1:num_gazes, 1:num_unique_points, 2:num_complete_repititions+1),4);

    %     % convert to stupid plotting convention for azimuth.
    %     % need this for gaussian functions or others that are not periodic on
    %     % the sphere.
    %     sel = unique_point_azimuth >=
    %     unique_point_azimuth(

    % convert data to be used for fitting to radians.
    unique_point_azimuth_r = unique_point_azimuth/180*pi;
    unique_point_elevation_r = unique_point_elevation/180*pi;
    unique_gaze_angle_r = unique_gaze_angle/180*pi;
    % in MOOG experiments, negative gaze angle is towards the positive z axis
    % and negative azimuth angle is counterclockwise (leftward) when viewed
    % from the positive z axis towards the origin.  This is opposite of the
    % spherical coordinates defined in Matlab (which we use for fitting).
    unique_gaze_angle_r = -unique_gaze_angle_r;

    if build_backdoor_load
        backdoor_load_file = fullfile(backdoor_dir, [FILE backdoor_load_ext]);
        save(backdoor_load_file, 'resp_mat', 'unique_stim_type', 'unique_gaze_angle_r', ...
            'unique_point_azimuth_r', 'unique_point_elevation_r', 'num_complete_repititions', ...
            'num_unique_points', 'fixation_type', 'num_gazes', 'unique_azimuth', ...
            'unique_elevation', 'unique_point_azimuth', 'unique_point_elevation', ...
            'unique_gaze_angle', 'max_stim_type');
        sprintf('data saved to %s', backdoor_load_file)
    end

else % ~use_backdoor_load
    backdoor_load_file = fullfile(backdoor_dir, [FILE backdoor_load_ext]);
    load(backdoor_load_file,'-mat');
end

% hack the data if we are only interested in fitting a single gaze angle
if fit_single_gaze_angle ~= CURVEFIT_INVALID_GAZE_ANGLE
    i = find(unique_gaze_angle == fit_single_gaze_angle);
    if isempty(i) | length(i) > 1
        error(['bad or multiple single gaze angle to fit in ' script_name]);
    end
    resp_mat = resp_mat(1:max_stim_type, i, 1:num_unique_points, 1:num_complete_repititions+1);
    unique_gaze_angle = unique_gaze_angle(i);
    unique_gaze_angle_r = unique_gaze_angle_r(i);
    num_gazes = 1;
end

% END load data

% START fit data
if ~run_fit, return; end
rand('state',sum(100*clock));   % seed it ups

% initialize global variable kludges that have to do with fitting.
curvefit_gaze_fits = zeros(num_gazes,num_unique_points);
curvefit_primary_fits = zeros(num_gazes,num_unique_points);
curvefit_secondary_fits = zeros(num_gazes,num_unique_points);
curvefit_gaze_residuals = zeros(num_gazes,num_unique_points);
curvefit_gaze_sse = zeros(1,num_gazes);
curvefit_gaze_mean_sse = zeros(1,num_gazes);

% iterate the curve fitting minimization some number of times.
% repeat for each model to be fit to.
num_stims = length(fitting_stims);
num_fits = length(fitting_models);
best_residual = -ones(num_stims,num_fits);
best_x = zeros(num_stims,num_fits,CURVEFIT_MAX_PARAMETERS);
all_residual = zeros(num_stims,num_fits,minimization_iterations);
num_params = zeros(num_stims,num_fits);
for i=1:num_stims

    % do not crash and burn if this data set does not contain the
    % stim type that we are interested in.
    if length( find(unique_stim_type == fitting_stims(i)) ) == 0
        sprintf('WARNING fit: %s does not contain stim %d data', ...
            FILE, fitting_stims(i))
        continue;
    end

    for j=1:num_fits

        if ~use_backdoor_fit  % there is no backdoor, you are not reading this.

            % reduce the dimension of the trial data for this stim type to 3.
            stim_data(1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1) = ...
                resp_mat(fitting_stims(i), 1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1);

            for n=1:minimization_iterations
                %sprintf('iter %d %s stim %d %s', n, func2str(fitting_models(j)),fitting_stims(i), FILE)

              matlab_version = version;
              if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                  % call the specified function to generate an initial guess
                  % for the parameters, and to create upper and lower bounds.
                  R = Curvefit_head_centered_bounds_half_ng(stim_data, num_gazes);
                  % run the curve fitting routine using the supplied constraints
                  [x,residual,exitflag,output] = fmincon( @Curvefit_head_centered_7p_asp_half_ng, R{CURVEFIT_BOUNDS_X0}, ...
                      [],[],[],[], R{CURVEFIT_BOUNDS_LB}, R{CURVEFIT_BOUNDS_UB},[], fitting_options, ...
                      @Curvefit_cos_tuning_7p_halfrectmodel, unique_point_azimuth_r, unique_point_elevation_r, ...
                      stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
                  sprintf('%d iterations %d funevals', output.iterations, output.funcCount)

              end
              
              if strcmp(matlab_version, '6.5.0.180913a (R13)')
                  % call the specified function to generate an initial guess
                  % for the parameters, and to create upper and lower bounds.
                  R = feval(fitting_bounds(j), stim_data, num_gazes);
                  
                  % run the curve fitting routine using the supplied constraints
                  [x,residual,exitflag,output] = fmincon( fitting_models(j), R{CURVEFIT_BOUNDS_X0}, ...
                      [],[],[],[], R{CURVEFIT_BOUNDS_LB}, R{CURVEFIT_BOUNDS_UB},[], fitting_options, ...
                      fitting_tuning_models(j), unique_point_azimuth_r, unique_point_elevation_r, ...
                      stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
                  sprintf('%d iterations %d funevals', output.iterations, output.funcCount)

              end
                
              
              
%                 Commented and then moved  underlying code to the to preceeding 'if' statements to accomoadate
%                 differences in matlab versions. -- Tunde
                
% %                 % call the specified function to generate an initial guess
% %                 % for the parameters, and to create upper and lower bounds.
% %                  R = feval(fitting_bounds(j), stim_data, num_gazes); % 
% % 
% %                 % run the curve fitting routine using the supplied constraints
% %                 [x,residual,exitflag,output] = fmincon( fitting_models(j), R{CURVEFIT_BOUNDS_X0}, ...
% %                     [],[],[],[], R{CURVEFIT_BOUNDS_LB}, R{CURVEFIT_BOUNDS_UB},[], fitting_options, ...
% %                     fitting_tuning_models(j), unique_point_azimuth_r, unique_point_elevation_r, ...
% %                     stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
% %                 sprintf('%d iterations %d funevals', output.iterations, output.funcCount)



                % keep around the best result over all iterations
                num_params(i,j) = length(R{CURVEFIT_BOUNDS_X0});
                if exitflag > 0
                    if best_residual(i,j) < 0 | residual < best_residual(i,j)
                        best_x(i,j,1:num_params(i,j)) = x;
                        best_residual(i,j) = residual;
                    end
                else
                    residual = -1;
                end
                all_residual(i,j,n) = residual;

                drawnow;  % I want to see!!!
            end

            if build_backdoor_fit
                fid=fopen(backdoor_fit_file, 'a');   % open output file
%                 commented section moved below and is now dependent on matlab version
%                 fprintf( fid, '%s %s %s %s %d %d %.16e ', FILE, func2str(fitting_models(j)), ...
%                     func2str(fitting_bounds(j)), func2str(fitting_tuning_models(j)), ...
%                     fitting_stims(i), num_params(i,j), best_residual(i,j) );
                
                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                fprintf( fid, '%s %s %s %s %d %d %.16e ', FILE, func2str(@Curvefit_head_centered_7p_asp_half_ng), ...
                    func2str(@Curvefit_head_centered_bounds_half_ng), func2str(@Curvefit_cos_tuning_7p_halfrectmodel), ...
                    fitting_stims(i), num_params(i,j), best_residual(i,j) );
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                fprintf( fid, '%s %s %s %s %d %d %.16e ', FILE, func2str(fitting_models(j)), ...
                    func2str(fitting_bounds(j)), func2str(fitting_tuning_models(j)), ...
                    fitting_stims(i), num_params(i,j), best_residual(i,j) );
                end
              
              
                fprintf( fid, '%.16e ', best_x(i,j,1:num_params(i,j)) );
                fprintf( fid, '| ' );
                fprintf( fid, '%e ', all_residual(i,j,:) );
                if ~plot_R2_export
                    fprintf( fid, '\n\n' );
                end
                fclose(fid);
            end

        else % ~use_backdoor_fit
            % read all the results from the backdoor fitting file
            [file_names fitting_model_str fitting_bound_str tuning_model_str ...
                file_stims file_num_params residual x0str all_residual_str] = ...
                textread( backdoor_fit_file, '%s %s %s %s %d %d %f %[^|]| %[^\n]' );
            num_files = length(file_names);

            % find the last result in the file that we are interested in
            select = ( strcmp(FILE, file_names) & ...
                strcmp(fitting_model_str, func2str(fitting_models(j))) & ...
                strcmp(fitting_bound_str, func2str(fitting_bounds(j))) & ...
                strcmp(tuning_model_str, func2str(fitting_tuning_models(j))) & ...
                file_stims == fitting_stims(i) );

            if ~any(select)
                str1 = sprintf( 'WARNING fit: backdoor file %s does not contain desired data for %s', ...
                    backdoor_fit_file, FILE );
                str2 = sprintf( '   %s %s %s stim %d', func2str(fitting_models(j)), ...
                    func2str(fitting_bounds(j)), func2str(fitting_tuning_models(j)), ...
                    fitting_stims(i) );
                [str1 str2]
                continue;
            else
                %                 % take the last result with the lowest residual
                %                 [m k] = sort(residual(select)); k = k(1);
                %                 select = diff([0; (cumsum(select) == k)]);
                %                 m = find(select == 1);
                %                 if length(m) ~= 1, error('something weird in ' script_name); end
                % take the last result in the file
                m = find(select == 1);
                m = m(length(m));
            end

            num_params(i,j) = file_num_params(m);
            best_residual(i,j) = residual(m);

            % extract the x0 parameters from the x0 string
            [best_x(i,j,1:num_params(i,j)), cnt] = sscanf( x0str{m}, '%f', inf );
            if cnt < num_params(i,j)
                ['Error reading x0 parameters "' x0str{m} '"']
                continue;
            end

            % extract all the residuals for each optimization iteration.
            % this causes matlab issues if different results have different
            % numbers of iterations, so read to a temp first before doing
            % the assign.
            [tmp, cnt] = sscanf( all_residual_str{m}, '%f', inf );
            if cnt < 1
                ['Error reading residuals "' all_residual_str{m} '"']
                continue;
            end
            all_residual(i,j,1:length(tmp)) = tmp;
        end

        %matrices for viewing seperated terms of plot ... hack for 1 fit for now
        %    [primary,secondary]=Curvefit_cos_tuning_7p_asp_seperator(best_x(1,1,1:num_params(i,j)),unique_point_azimuth, unique_point_elevation);
    end % for each fit
end % for each stim type

% END fit data

% START compute statistics
if ~run_compute_stat, return; end

%compute d'


% compute values from the data and the best fit data:
SS_YY = zeros(num_stims);
SSE = zeros(num_stims,num_fits);
R_2 = zeros(num_stims,num_fits);
bootstrap_R_2_distrib = zeros(num_stims,num_fits,bootstrap_num_samples);
% in order of dimensions -> stim, fit, gaze, point
% NOTE: this is different than the format of the raw data,
% but it is confusing since they are both for dimensional.
% stim type is indexed by unique_stim_type index, and not by unique stim type.
fit_data = zeros(num_stims,num_fits,num_gazes,num_unique_points);
fit_residual = zeros(num_stims,num_fits,num_gazes,num_unique_points);
primary_data = zeros(num_stims,num_fits,num_gazes,num_unique_points);
secondary_data = zeros(num_stims,num_fits,num_gazes,num_unique_points);
for i=1:num_stims

    % do not crash and burn if this data set does not contain the
    % stim type that we are interested in.
    if length( find(unique_stim_type == fitting_stims(i)) ) == 0
        sprintf('WARNING stat: %s does not contain stim %d data', ...
            FILE, fitting_stims(i))
        continue;
    end

    % reduce the dimension of the trial data for this stim type to 3.
    stim_data(1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1) = ...
        resp_mat(fitting_stims(i), 1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1);
    mean_trial_data(1:num_gazes, 1:num_unique_points) = ...
        resp_mat(fitting_stims(i), 1:num_gazes, 1:num_unique_points, 1);
    trial_data(1:num_gazes, 1:num_unique_points, 2:num_complete_repititions+1) = ...
        resp_mat(fitting_stims(i),1:num_gazes, 1:num_unique_points, 2:num_complete_repititions+1);

    % compute the mean over all repititions, gaze angles, and points.
    % do this for the sqrt of the data, since that is how we did the fitting.
    mean_data = sum(sum(sum(Curvefit_sqrt_abs(trial_data)))) / ...
        (num_gazes * num_complete_repititions * num_unique_points);

    % compute the ss_yy term of the correlation coefficient.
    % this is the variance of the data times N
    % (in our case N = num_gazes * num_complete_repititions).
    % sum variance over all gazes, repitions, and points to get total
    % ss_yy
    % again, do this for the sqrt of the data, since we fit that way.
    SS_YY(i) = sum(sum(sum((Curvefit_sqrt_abs(trial_data) - mean_data).^2)));

    % compute ss_yy using the mean data at each point, instead of trial data
    mean_data = sum(sum(Curvefit_sqrt_abs(mean_trial_data))) / ...
        (num_gazes * num_unique_points);
    SS_YY(i) = sum(sum((Curvefit_sqrt_abs(mean_trial_data) - mean_data).^2));

    for j=1:num_fits

        % do not crash and burn if we do not have a best fit.
        if best_residual(i,j) == -1
            sprintf('WARNING stat: %s does not contain stim %d fit %s data', ...
                FILE, fitting_stims(i), func2str(fitting_models(j)))
            continue;
        end

% %         moved this section to depend on the matlab version
%         % compute SSE and verify with fitting results residual.
%         SSE(i,j) = feval( fitting_models(j), best_x(i,j,1:num_params(i,j)), ...
%             fitting_tuning_models(j), unique_point_azimuth_r, unique_point_elevation_r, ...
%             stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
        
        
        if strcmp(matlab_version, '7.2.0.232 (R2006a)')
            % compute SSE and verify with fitting results residual.
            SSE(i,j) = Curvefit_head_centered_7p_asp_half_ng( best_x(i,j,1:num_params(i,j)), ...
                @Curvefit_cos_tuning_7p_halfrectmodel, unique_point_azimuth_r, unique_point_elevation_r, ...
                stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
        end

        if strcmp(matlab_version, '6.5.0.180913a (R13)')
            % compute SSE and verify with fitting results residual.
            SSE(i,j) = feval( fitting_models(j), best_x(i,j,1:num_params(i,j)), ...
                fitting_tuning_models(j), unique_point_azimuth_r, unique_point_elevation_r, ...
                stim_data, num_complete_repititions, unique_gaze_angle_r, fixation_type );
        end

        e = 0.5;
        if SSE(i,j) < best_residual(i,j) - e | SSE(i,j) > best_residual(i,j) + e
            %error( ['can not reproduce SSE in ' script_name] );
        end

        % save fits at each gaze angle.
        fit_data(i,j,:,:) = curvefit_gaze_fits;
        % save residuals versus the mean at each gaze angle.
        fit_residual(i,j,:,:) = curvefit_gaze_residuals;
        % save sse over all trials at each gaze angle.
        fit_sse(i,j,:) = curvefit_gaze_sse;
        primary_data(i,j,:,:) = curvefit_primary_fits;
        secondary_data(i,j,:,:) = curvefit_secondary_fits;

        % compute SSE using the mean data at each point, instead of trial data
        SSE(i,j) = sum(curvefit_gaze_mean_sse);
        fit_sse(i,j,:) = curvefit_gaze_mean_sse;

        % compute the r^2 correlation coefficient.
        R_2(i,j) = 1 - SSE(i,j) / SS_YY(i);

        % to bootstrap?
        if compute_stat_bootstrap & num_complete_repititions < 2
            sprintf('WARNING stat: %s only contains %d reptitions, skipping bootstrap', ...
                FILE, num_complete_repititions)
        elseif compute_stat_bootstrap
            % roll up fit data into a row vector with points for a single
            % gaze angle next to each other.
            curvefit_bootstrap_fits = curvefit_gaze_fits';
            curvefit_bootstrap_fits = curvefit_bootstrap_fits(:)';
            % roll up trial data in the same fashion so that the callback
            % function does not have to rearrange anything.  bootstrap
            % requires that the data to be resampled be along the rows.
            bootstrap_data = permute(trial_data,[3,2,1]);
            bootstrap_data = bootstrap_data(:,:);
            bootstrap_R_2_distrib(i,j,:) = bootstrp(bootstrap_num_samples, ...
                'Curvefit_bootstrap_R_2', bootstrap_data);
            %sprintf( '%s %.10e', func2str(fitting_models(j)), ...
            %    max(bootstrap_R_2_distrib(i,j,:)) )
        end
    end % for each fit
end % for each stim type

%compute chi2 for stim 2 only- Z. Briggs
%convert data to columns for chi2
%**** GCD: this section only works properly right now if num_stims = 1
chix=[];
for j=1:num_complete_repititions
    for i=1:num_unique_points
        chix=[chix;i];
    end
end

chiy=[];
for i=2:num_complete_repititions+1  %first value is the mean
    for j=1:num_unique_points
        chiy=[chiy;resp_mat(fitting_stims(1),1,j,i)];
    end
end
chi2_stat=zeros(num_fits,1);
chi2p_stat=zeros(num_fits,1);
for c=1:num_fits
    [chi2_stat(c),chi2p_stat(c)] = Chi2_Test_Fval(chix, chiy, fit_data(1,c,1,:), num_params(c));
end;

% END compute statistics

% START plot data
if ~run_plot, return; end
if plot_contour

    % convert the azimuth and elevation into a format amenable to plotting
    unique_azimuth = sort(unique_azimuth)' - 90;  % is the sort redundant?
    unique_elevation = sort(unique_elevation)';
    unique_azimuth = [unique_azimuth unique_azimuth(1) + 360];  % wrap around sphere
    [plot_azimuth plot_elevation] = meshgrid(unique_azimuth, unique_elevation);
    num_grid_azimuth = length(unique_azimuth);
    num_grid_elevation = length(unique_elevation);
    plot_mean_data = zeros(num_grid_elevation,num_grid_azimuth);
    plot_fit_data = zeros(num_grid_elevation,num_grid_azimuth);
    plot_residual_data = zeros(num_grid_elevation,num_grid_azimuth);
    plot_primary_data = zeros(num_grid_elevation,num_grid_azimuth);
    plot_secondary_data = zeros(num_grid_elevation,num_grid_azimuth);

    % set the bounds and ticks for the axis on the contour plots,
    % depending on how we formatted the data for plotting above.
    x_min = -90; x_max = 270; y_min = -90; y_max = 90;
    x_tick = [-90:45:270]; y_tick = [-90:45:90];

    for i=1:num_stims

        % do not crash and burn if this data set does not contain the
        % stim type that we are interested in.
        if length( find(unique_stim_type == fitting_stims(i)) ) == 0
            sprintf('WARNING plot: %s does not contain stim %d data', ...
                FILE, fitting_stims(i))
            continue;
        end

        % reduce dim for the data for this stim type
        stim_data(1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1) = ...
            resp_mat(fitting_stims(i), 1:num_gazes, 1:num_unique_points, 1:num_complete_repititions+1);
        % this is a kludge to make sure singular dimensions are not collapsed
        stim_fit_data = zeros(10,10,num_unique_points);
        stim_primary_data = zeros(10,10,num_unique_points);
        stim_secondary_data = zeros(10,10,num_unique_points);
        stim_residual = zeros(10,10,num_unique_points);
        m = size(fit_data,2); n = size(fit_data,3);
        stim_fit_data(1:m,1:n,:) = fit_data(i,:,:,:);
        stim_primary_data(1:m,1:n,:) = primary_data(i,:,:,:);
        stim_secondary_data(1:m,1:n,:) = secondary_data(i,:,:,:);
        stim_residual(1:m,1:n,:) = fit_residual(i,:,:,:);

        % get the min, max, and range of the data and the fitted data
        % over repititions, gaze angles, and points
        max_data = max( [max(max(max(stim_data))) max(max(max(stim_fit_data)))] );
        min_data = min( [min(min(min(stim_data))) min(min(min(stim_fit_data)))] );
        range_data = max_data - min_data;

        % specify the elevations to draw contour lines at for data and residual
        % plots.  specifying the number of contour lines to draw does not take
        % the min and max range of all the data into account, so it is
        % difficult to compare between contours.
        data_clines = 0.99*[min_data:range_data/10:max_data]  % 0.99 to see max and min points
        residual_clines = 0.99*[-range_data/2:range_data/10:range_data/2]

        % this is something of a hack and mostly intended for three fits.
        %for j=1:num_fits
        for j=1:plot_contour_compares_step:plot_contour_num_compares

            % do not crash and burn if we do not have a best fit.
            if best_residual(i,j) == -1
%                 sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
%                     FILE, fitting_stims(i), func2str(fitting_models(j)))
                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
                        FILE, fitting_stims(i), func2str(@Curvefit_head_centered_7p_asp_half_ng))
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
                        FILE, fitting_stims(i), func2str(fitting_models(j)))
                end
                                
                continue;
            end

            % make plots to compare each fit with the following fit, and the last
            % fit with the first fit.
            k = mod(j,num_fits) + 1;

            % create a figure
%                    hmain = figure(num_fits*(i-1) + j + start_figure_contour - 1); h = hmain;
%                   clf reset;  % clear the figure
            ss = get(0,'ScreenSize');
            pos = [0.2*ss(3) 0.1*ss(4) 0.8*ss(3) 0.8*ss(4)];
            %set(h,'Position',pos);

            % string to identify direction of gaze angle change
            if fixation_type == CURVEFIT_VARY_FIXATION_X
                fixation_type_str = 'x';
            elseif fixation_type == CURVEFIT_VARY_FIXATION_Y
                fixation_type_str = 'y';
            end

            % these plots are hard coded to only have three gaze angles and two fits.
            % there is a semi-kludge for plotting one gaze angle with one fit.
            pstr = {'' '' ''};
            max_3_gazes = min([3 num_gazes]);
            max_2_fits = min([3 num_fits+1]);
            for gas=1:max_3_gazes

                % convert the mean and fit data into a format amenable to plotting
                for m=1:num_grid_elevation
                    for n=1:num_grid_azimuth
                        % i can not think of a cleaner way to do this stupid
                        % all azimuths at the pole thing.
                        e = 1e-10;  % floating point tolerance
                        select = ((abs(unique_point_azimuth - mod(unique_azimuth(n),360)) < e | ...
                            unique_elevation(m) == 90 | unique_elevation(m) == -90) & ...
                            abs(unique_point_elevation - unique_elevation(m)) < e);
                        if sum(select) == 0;
                            eat_shit_and_die = 1;
                        end
                        plot_mean_data(m,n) = stim_data(gas,select,1);
                        plot_fit_data(m,n) = stim_fit_data(j,gas,select);
                        plot_primary_data(m,n) = stim_primary_data(j,gas,select);
                        plot_secondary_data(m,n) = stim_secondary_data(j,gas,select);
                        plot_residual_data(m,n) = stim_residual(j,gas,select);
                        if num_fits > 1
                            plot_fit_cmp_data(m,n) = stim_fit_data(k,gas,select);
                            plot_residual_cmp_data(m,n) = stim_residual(k,gas,select);
                            plot_primary_cmp_data(m,n) = stim_primary_data(k,gas,select);
                            plot_secondary_cmp_data(m,n) = stim_secondary_data(k,gas,select);
                        end
                    end
                end


                %New Sphere Plot
                %                 plot_sphere_data=zeros(26);
                %                 for m=1:num_unique_points
                %                     for n=1:num_unique_points
                %                     plot_sphere_data(m,n)=Curvefit_cos_tuning_7p_gauss(best_x,unique_point_azimuth_r(n),unique_point_elevation_r(m));
                %                     sphplot(m,n)=sph2cart(unique_point_azimuth_r(n),unique_point_elevation(m),plot_sphere_data(m,n));
                %                     end
                %                 end
                %
                %                 plot3(unique_point_azimuth_r,unique_point_elevation_r,plot_sphere_data)
                %


                %plot primary and secondary data
                z=figure;
                set(gcf,'Position',pos);

                colormap('gray');
                cmap = colormap;
                colormap(1-cmap);
                subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+3);
                %subplot('Position',[0 0.5 0.5 0.5]);
                [C h] = contourf(plot_azimuth,sin(plot_elevation*pi/180),plot_primary_data);
                colorbar;
                set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                %axis([x_min,x_max,y_min,y_max]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                axis([x_min,x_max,-1, 1]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                %caxis([-7, 7]);
                %xlabel('azimuth'); ylabel('elevation');
                
% %                 below code now runs depending on version
%                 str = {sprintf('%s', func2str(fitting_models(j)))};                                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = {sprintf('%s', func2str(fitting_models(j)))};
                end
                             
                if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end

                               
                subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+4);
                %subplot('Position',[0 0.5 0.5 0.5]);
                [C h] = contourf(plot_azimuth,sin(plot_elevation*pi/180),plot_secondary_data);

                set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                %axis([x_min,x_max,y_min,y_max]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                axis([x_min,x_max,-1,1]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                %caxis([-7, 7]);
                %xlabel('azimuth'); ylabel('elevation');
                colorbar;
%                 str = {sprintf('%s \n Weighting 180 = %d', func2str(fitting_models(j)), best_x(1,1,7))};
                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = {sprintf('%s \n Weighting 180 = %d', func2str(@Curvefit_head_centered_7p_asp_half_ng), best_x(1,1,7))};
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = {sprintf('%s \n Weighting 180 = %d', func2str(fitting_models(j)), best_x(1,1,7))};
                end
                
                
                if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end


                if num_fits>1
                    subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+1);
                    %subplot('Position',[0 0.5 0.5 0.5]);
                    [C h] = contourf(plot_azimuth,plot_elevation,plot_primary_cmp_data);
                    colorbar;
                    set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                    axis([x_min,x_max,y_min,y_max]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                    %xlabel('azimuth'); ylabel('elevation');
                    
%                     str = {sprintf('%s', func2str(fitting_models(k)))};
                    if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                        str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                    end

                    if strcmp(matlab_version, '6.5.0.180913a (R13)')
                        str = {sprintf('%s', func2str(fitting_models(k)))};
                    end
                    
                    if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end

                    subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+2);
                    %subplot('Position',[0 0.5 0.5 0.5]);
                    [C h] = contourf(plot_azimuth,plot_elevation,plot_secondary_cmp_data);
                    colorbar;
                    set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                    axis([x_min,x_max,y_min,y_max]); caxis('auto'); if plot_clabels, clabel(C,h), end;
                    %xlabel('azimuth'); ylabel('elevation');

%                     str = {sprintf('%s \n Weighting 180 = %d', func2str(fitting_models(k)), best_x(1,2,6))};                        
                    if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                        str = {sprintf('%s \n Weighting 180 = %d', func2str(@Curvefit_head_centered_7p_asp_half_ng), best_x(1,2,6))};
                    end

                    if strcmp(matlab_version, '6.5.0.180913a (R13)')
                        str = {sprintf('%s \n Weighting 180 = %d', func2str(fitting_models(k)), best_x(1,2,6))};
                    end
                    
                    if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end
                end






                pause on;
                y=figure;
                set(gcf,'Position',pos);

        % get the min, max, and range of the data and the fitted data
        % mean responses
        max_mean_data = max( [max(max(plot_mean_data)) max(max(plot_fit_data))] );
        min_mean_data = min( [min(min(plot_mean_data)) min(min(plot_fit_data))] );
        range_mean_data = max_mean_data - min_mean_data;
        r_clines = [-range_mean_data/2:range_mean_data/10:range_mean_data/2]
        colormap('gray');
        cmap = colormap;
        colormap(1-cmap);
        
                % plot the mean data
                subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+1);
                %subplot('Position',[0 0 0.5 0.5]);
                [C h] = contourf(plot_azimuth,sin(plot_elevation*pi/180),plot_mean_data,data_clines);
                colorbar;
                set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                %axis([x_min,x_max,y_min,y_max]); caxis([min_mean_data max_mean_data]); if plot_clabels, clabel(C,h), end;
                axis([x_min,x_max,-1,1]); caxis([min_mean_data max_mean_data]); if plot_clabels, clabel(C,h), end;
                xlabel('azimuth'); ylabel('elevation');
                str = {sprintf('%s', FILE) sprintf('gaze_angle=%d vary_fixation_%s', ...
                    unique_gaze_angle(gas), fixation_type_str)};
                h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth');

                % plot the first fit
                subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+2);
                %subplot('Position',[0 0.5 0.5 0.5]);
                [C h] = contourf(plot_azimuth,sin(plot_elevation*pi/180),plot_fit_data,data_clines);
                colorbar;
                set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                %axis([x_min,x_max,y_min,y_max]); caxis([min_mean_data max_mean_data]); if plot_clabels, clabel(C,h), end;
                axis([x_min,x_max,-1,1]); caxis([min_mean_data max_mean_data]); if plot_clabels, clabel(C,h), end;
                %xlabel('azimuth'); ylabel('elevation');
                
%                 str = {sprintf('%s', func2str(fitting_models(j)))};                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = {sprintf('%s', func2str(fitting_models(j)))};
                end
                                                    
                if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end

                % plot the first residual
                subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+max_2_fits+2);
                %subplot('Position',[0 0 0.5 0.5]);
                %[C h] = contourf(plot_azimuth,plot_elevation,plot_residual_data, [-.5 .5 -2 2 -4 4]);
                [C h] = contourf(plot_azimuth,sin(plot_elevation*pi/180),(plot_residual_data+0.5*range_mean_data), data_clines);
                colorbar;
                set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                %axis([x_min,x_max,y_min,y_max]); caxis([min_mean_data max_mean_data]); 
                axis([x_min,x_max,-1,1]); caxis([min_mean_data max_mean_data]); 
                if plot_clabels, clabel(C,h), end
                %xlabel('azimuth'); ylabel('elevation');
                
%                 str = [sprintf('%s SSE=%.5f ', func2str(fitting_models(j)), fit_sse(i,j,gas)) ];                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = [sprintf('%s SSE=%.5f ', func2str(@Curvefit_head_centered_7p_asp_half_ng), fit_sse(i,j,gas)) ];
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = [sprintf('%s SSE=%.5f ', func2str(fitting_models(j)), fit_sse(i,j,gas)) ];
                end                
   
                str = [str repmat([' '], 1, plot_str_len - size(str,2))];
                pstr{gas} = [pstr{gas};str];

                if num_fits > 1
                    % plot the second fit
                    subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+3);
                    %subplot('Position',[0 0.5 0.5 0.5]);
                    [C h] = contourf(plot_azimuth,plot_elevation,plot_fit_cmp_data,data_clines);
                    set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                    axis([x_min,x_max,y_min,y_max]); caxis([min_data max_data]); if plot_clabels, clabel(C,h), end;
                    %xlabel('azimuth'); ylabel('elevation');
                    
%                     str = {sprintf('%s', func2str(fitting_models(k)))};
                    if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                        str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                    end

                    if strcmp(matlab_version, '6.5.0.180913a (R13)')
                        str = {sprintf('%s', func2str(fitting_models(k)))};
                    end                  
                    
                    if gas == 1, h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); end

                    % plot the second residual
                    subplot(2*max_3_gazes,max_2_fits,2*max_2_fits*(gas-1)+6);
                    %subplot('Position',[0 0 0.5 0.5]);
                    [C h] = contourf(plot_azimuth,plot_elevation,plot_residual_cmp_data,residual_clines);
                    set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                    axis([x_min,x_max,y_min,y_max]); caxis([-range_data/2 range_data/2]); if plot_clabels, clabel(C,h), end;
                    %xlabel('azimuth'); ylabel('elevation');
                    
                    
%                     str = [sprintf('%s SSE=%.5f ', func2str(fitting_models(k)), fit_sse(i,k,gas)) ];                    
                    if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                        str = [sprintf('%s SSE=%.5f ', func2str(@Curvefit_head_centered_7p_asp_half_ng), fit_sse(i,k,gas)) ];
                    end

                    if strcmp(matlab_version, '6.5.0.180913a (R13)')
                        str = [sprintf('%s SSE=%.5f ', func2str(fitting_models(k)), fit_sse(i,k,gas)) ];
                    end

                    str = [str repmat([' '], 1, plot_str_len - size(str,2))];
                    pstr{gas} = [pstr{gas};str];
                end

                %pstr{gas} = {pstr{gas} 'happy day'};
            end

            % string formatting sucks.  you figure it out.
            str = sprintf('%-70s', ' ');  % len must be the same as plot_str_len
            for gas=max_3_gazes+1:3
                pstr{gas} = [pstr{gas};str;str];
            end
            for gas=1:2
                pstr{gas} = [pstr{gas};str;str;str];
            end
            pstr{max_3_gazes} = [pstr{max_3_gazes};str];
            
%             str = sprintf( '%s R_2=%-.5f', func2str(fitting_models(j)), R_2(i,j) );     
            if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                str = sprintf( '%s R_2=%-.5f', func2str(@Curvefit_head_centered_7p_asp_half_ng), R_2(i,j) );
            end

            if strcmp(matlab_version, '6.5.0.180913a (R13)')
                str = sprintf( '%s R_2=%-.5f', func2str(fitting_models(j)), R_2(i,j) );
            end
      
            str = [str repmat([' '], 1, plot_str_len - size(str,2))];
            pstr{max_3_gazes} = [pstr{max_3_gazes};str];
            if num_fits > 1
%                 str = sprintf( '%s R_2=%-.5f', func2str(fitting_models(k)), R_2(i,k) );              
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = sprintf( '%s R_2=%-.5f', func2str(@Curvefit_head_centered_7p_asp_half_ng), R_2(i,k) );
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = sprintf( '%s R_2=%-.5f', func2str(fitting_models(k)), R_2(i,k) );
                end                               
                
                str = [str repmat([' '], 1, plot_str_len - size(str,2))];
                pstr{max_3_gazes} = [pstr{max_3_gazes};str];
            end

            %      for gas=1:3
            %         if max_3_gazes == 3
            %             axes('Position',[.02 0.71-0.285*(gas-1) 1 1]);
            %         elseif max_3_gazes == 2
            %             axes('Position',[.02 0.67-0.47*(gas-1) 1 1]);
            %         else
            %             axes('Position',[.02 0.4-0.285*(gas-1) 1 1]);
            %         end
            gas=1
            axes('Position',[.02 0.4-0.285*(gas-1) 1 1]);     %was .285(gas-1)
            h = text(0,0,pstr{gas});
            set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth','FontSize',10);



            %plot chi stats
            for c=1:num_fits
%                 str=sprintf('%s \nChi2 = %f \nChi2P = %f \nnum_complete_reps = %f', func2str(fitting_models(c)),chi2_stat(c), chi2p_stat(c), num_complete_repititions);            
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                str=sprintf('%s \nChi2 = %f \nChi2P = %f \nnum_complete_reps = %f', func2str(@Curvefit_head_centered_7p_asp_half_ng),chi2_stat(c), chi2p_stat(c), num_complete_repititions);
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                str=sprintf('%s \nChi2 = %f \nChi2P = %f \nnum_complete_reps = %f', func2str(fitting_models(c)),chi2_stat(c), chi2p_stat(c), num_complete_repititions);
                end
                
                h = text(0,-.045-.05*(k+c),str);
                set(h, 'Interpreter', 'none'); set(h,'FontName','FixedWidth','FontSize',10);
            end

            axis off;

            % optionally save the figure in the backdoor directory
            if save_contour_figures

%                 figure_name = fullfile(backdoor_dir, [FILE '_' func2str(fitting_models(j)) '_to_' ...
%                     func2str(fitting_models(k)) sprintf('_stim%d', i) save_figure_ext]);                                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    figure_name = fullfile(backdoor_dir, [FILE '_' func2str(@Curvefit_head_centered_7p_asp_half_ng) '_to_' ...
                        func2str(@Curvefit_head_centered_7p_asp_half_ng) sprintf('_stim%d', i) save_figure_ext]);
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    figure_name = fullfile(backdoor_dir, [FILE '_' func2str(fitting_models(j)) '_to_' ...
                        func2str(fitting_models(k)) sprintf('_stim%d', i) save_figure_ext]);
                end                     

%                 saveas(hmain, figure_name, save_figure_format); %Originally -- BA edit
                saveas(h, figure_name, save_figure_format);
            end

            % optionally print a hard copy
            if print_contour_figures
                print(printer_contour_figures);
            end

        end % for each fit
    end % for each stim type

end % if plot_contour

if plot_R2_export
    curvefit_R2_export = R_2;
    curvefit_Chi2_export = chi2p_stat(1);
    curvefit_fitting_stims_export = fitting_stims;
    if build_backdoor_fit
        fid=fopen(backdoor_fit_file, 'a');   % open output file
        fprintf( fid ,'| ');
        fprintf( fid ,'%f %f ',curvefit_R2_export, curvefit_Chi2_export);
        fprintf( fid , '\n\n');
        fclose(fid);
    end

    for i=1:num_stims
        if length( find(unique_stim_type == fitting_stims(i)) ) == 0
            curvefit_fitting_stims_export(i) = -1;
        end
    end
end

if plot_param_export
    curvefit_gaze_angle_export = unique_gaze_angle;
    curvefit_fixation_type_export = fixation_type;
    curvefit_fitting_param_export = best_x;
    curvefit_fitting_stims_export = fitting_stims;
    for i=1:num_stims
        if length( find(unique_stim_type == fitting_stims(i)) ) == 0
            curvefit_fitting_stims_export(i) = -1;
        end
    end
end

% see Curvefit_defines for enumeration types.
curvefit_R2_distrib_sig_export = ...
    CURVEFIT_R2_DISTRIB_EQUAL_MEANS*ones(num_stims, ...
    ceil(plot_contour_num_compares/plot_contour_compares_step));

if plot_R2_distrib_export & num_complete_repititions < 2
    sprintf('WARNING plot: %s only contains %d reptitions, skipping bootstrap', ...
        FILE, num_complete_repititions)
elseif plot_R2_distrib_export

    % sort all the r^2 distributions so we can get the confidence intervals below.
    bootstrap_R_2_distrib_sorted = sort(bootstrap_R_2_distrib,3);

    for i=1:num_stims

        % do not crash and burn if this data set does not contain the
        % stim type that we are interested in.
        if length( find(unique_stim_type == fitting_stims(i)) ) == 0
            sprintf('WARNING plot: %s does not contain stim %d data', ...
                FILE, fitting_stims(i))
            continue;
        end

        % this is something of a hack and mostly intended for three fits.
        %for j=1:num_fits
        for j=1:plot_contour_compares_step:plot_contour_num_compares

            % do not crash and burn if we do not have a best fit.
            if best_residual(i,j) == -1
%                 sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
%                     FILE, fitting_stims(i), func2str(fitting_models(j)))              
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
                        FILE, fitting_stims(i), func2str(@Curvefit_head_centered_7p_asp_half_ng))
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    sprintf('WARNING plot: %s does not contain stim %d fit %s data', ...
                        FILE, fitting_stims(i), func2str(fitting_models(j)))
                end      
                
                continue;
            end

            % make plots to compare each fit with the following fit, and the last
            % fit with the first fit.
            k = mod(j,num_fits) + 1;

            % decide if the better fit is statistically significant by
            % checking if the 5% confidence from the r^2 distrib with the
            % greater mean is greater than the 95% confidence from the r^2
            % distrib with the lesser mean.  this is mildly ugly code, but
            % i'm feeling lazy.  please make it cleaner for me.
            R_2_conf_j(1) = bootstrap_R_2_distrib_sorted(i,j,0.05*bootstrap_num_samples);
            R_2_conf_j(2) = bootstrap_R_2_distrib_sorted(i,j,0.95*bootstrap_num_samples);
            R_2_conf_k(1) = bootstrap_R_2_distrib_sorted(i,k,0.05*bootstrap_num_samples);
            R_2_conf_k(2) = bootstrap_R_2_distrib_sorted(i,k,0.95*bootstrap_num_samples);
            tmp = bootstrap_R_2_distrib(i,j,:);
            r2std_j = std(tmp); r2mean_j = mean(tmp);
            tmp = bootstrap_R_2_distrib(i,k,:);
            r2std_k = std(tmp); r2mean_k = mean(tmp);
            if r2mean_j > r2mean_k
                
%                 r2_sig_str = func2str(fitting_models(j));                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    r2_sig_str = func2str(@Curvefit_head_centered_7p_asp_half_ng);
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    r2_sig_str = func2str(fitting_models(j));
                end

                curvefit_R2_distrib_sig_export(i,j) = CURVEFIT_R2_DISTRIB_1_NSIG;
                if R_2_conf_j(1) >= R_2_conf_k(2)
                    r2_sig_str = [r2_sig_str ' YES'];
                    curvefit_R2_distrib_sig_export(i,j) = CURVEFIT_R2_DISTRIB_1_SIG;
                else
                    r2_sig_str = [r2_sig_str ' NO'];
                end
            elseif r2mean_k > r2mean_j
%                 r2_sig_str = func2str(fitting_models(k));        
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    r2_sig_str = func2str(@Curvefit_head_centered_7p_asp_half_ng);
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    r2_sig_str = func2str(fitting_models(k));
                end
                                
                curvefit_R2_distrib_sig_export(i,j) = CURVEFIT_R2_DISTRIB_2_NSIG;
                if R_2_conf_k(1) >= R_2_conf_j(2)
                    r2_sig_str = [r2_sig_str ' YES'];
                    curvefit_R2_distrib_sig_export(i,j) = CURVEFIT_R2_DISTRIB_2_SIG;
                else
                    r2_sig_str = [r2_sig_str ' NO'];
                end
            else
                r2_sig_str = 'EQUAL means';
                curvefit_R2_distrib_sig_export(i,j) = CURVEFIT_R2_DISTRIB_EQUAL_MEANS;
            end

            if plot_R2_distrib
                % create a figure
%                 hmain = figure(num_fits*(i-1) + j + start_figure_R2_distrib - 1); h = hmain;
%                 clf reset;  % clear the figure
                %             ss = get(0,'ScreenSize');
                %             pos = [0.2*ss(3) 0.1*ss(4) 0.8*ss(3) 0.8*ss(4)];
                %             set(h,'Position',pos);

                % a pseudo-option, 4 standard deviations makes a good looking hist.
                % changed this to use conf intervals as limits and for tick marks.
                %nstd = 4;
                %set(gca,'xtick',[r2mean-nstd*r2std:nstd/2*r2std:r2mean+nstd*r2std]);

                subplot(1,2,1);
                tmp = bootstrap_R_2_distrib_sorted(i,j,:);
                hist(tmp, plot_R2_distrib_num_bins);
                x_tick = [R_2_conf_j(1) r2mean_j R_2_conf_j(2)];
                set(gca,'xtick',x_tick,'xticklabel',sprintf('%.3f|',x_tick));
                xlim([r2mean_j-4*r2std_j r2mean_j+2*r2std_j]);
%                 str = {sprintf('%s', func2str(fitting_models(j)))};                                
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = {sprintf('%s', func2str(fitting_models(j)))};
                end                
                
                h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth');
                ylabel(FILE);
                h = xlabel(r2_sig_str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth');

                subplot(1,2,2);
                tmp = bootstrap_R_2_distrib_sorted(i,k,:);
                hist(tmp, plot_R2_distrib_num_bins);
                x_tick = [R_2_conf_k(1) r2mean_k R_2_conf_k(2)];
                set(gca,'xtick',x_tick,'xticklabel',sprintf('%.3f|',x_tick));
                xlim([r2mean_k-4*r2std_k r2mean_k+2*r2std_k]);
%                 str = {sprintf('%s', func2str(fitting_models(k)))};       
                if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                    str = {sprintf('%s', func2str(@Curvefit_head_centered_7p_asp_half_ng))};
                end

                if strcmp(matlab_version, '6.5.0.180913a (R13)')
                    str = {sprintf('%s', func2str(fitting_models(k)))};
                end
                                
                h = title(str); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth');

                % optionally save the figure in the backdoor directory
                if save_R2_distrib_figures
%                     figure_name = fullfile(backdoor_dir, [FILE '_r2_' func2str(fitting_models(j)) '_' ...
%                         func2str(fitting_models(k)) sprintf('_stim%d', i) save_figure_ext]);                   
                    if strcmp(matlab_version, '7.2.0.232 (R2006a)')
                        figure_name = fullfile(backdoor_dir, [FILE '_r2_' func2str(@Curvefit_head_centered_7p_asp_half_ng) '_' ...
                            func2str(@Curvefit_head_centered_7p_asp_half_ng) sprintf('_stim%d', i) save_figure_ext]);
                    end

                    if strcmp(matlab_version, '6.5.0.180913a (R13)')
                        figure_name = fullfile(backdoor_dir, [FILE '_r2_' func2str(fitting_models(j)) '_' ...
                            func2str(fitting_models(k)) sprintf('_stim%d', i) save_figure_ext]);
                    end
                                      
                    saveas(hmain, figure_name, save_figure_format);
                end

                % optionally print a hard copy
                if print_R2_distrib_figures
                    print(printer_R2_distrib_figures);
                end
            end % if plot R2 distrib

        end % for each plot
    end % for each stim type
end % if plot_R2_distrib
% END plot datapause;
% pause on;
% pause;
% close(z);
% close(y);
% pause off;

return;

%-----------------------------------------------------------------------------------------------------------------------
% END OF LINE
%-----------------------------------------------------------------------------------------------------------------------
