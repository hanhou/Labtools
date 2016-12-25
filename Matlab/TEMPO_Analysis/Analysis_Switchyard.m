%----------------------------------------------------------------------------------------------------------
%-- AnalysisSwitchyard.m: This function takes in the data, protocol, etc. and chains off to the
%-- 	relevant protocol-specific analysis code.  GCD, 1/3/2000
%--     Last Revised BJP, 3/1/02
%----------------------------------------------------------------------------------------------------------
function Analysis_Switchyard(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetTime, StopOffsetTime, PATH, FILE, batch_flag, UseSyncPulses);

global JustEditIt;

TEMPO_Defs;
% ProtocolDefs;
Path_Defs;	%this sets up paths to the analysis funtions, both Common and ProtocolSpecific

if(~isempty(batch_flag))
    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    t=sprintf('Doing analysis (%s) for Protocol: %s', Analysis{1}, protocol_names{Protocol+1});
    set(ListHandle, 'String', t);
end

%following function checks offset times and outputs starting and ending analysis (in spike bins)
num_trials = size(data.event_data, 3);

if EndTrial == -1
    EndTrial = num_trials;
end

[StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffsetTime, StopOffsetTime, UseSyncPulses);
data.UseSyncPulses=UseSyncPulses;


if (~isempty(data.spike_data))
    %compute the firing rate over all trials during the period between StartCode and StopCode
    data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
    
    % Plot firing rate against trials that will be used for further analysis. HH20150407
    if(isempty(batch_flag))
        
        % If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
        % Else, if all elements are negative, they are trials to be excluded.
        % This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
        select_trials = false(num_trials,1);
        if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
            select_trials(BegTrial:EndTrial) = true;
        elseif all(BegTrial > 0) % To be included
            select_trials(BegTrial) = true;
        elseif all(BegTrial < 0) % To be excluded
            select_trials(-BegTrial) = true;
            select_trials = ~ select_trials;
        else
            disp('Trial selection error...');
            keyboard;
        end
        
        % Plot firing rate as a function of trial number (from START to END of each trial) to check recording stability
        if ~isempty(SpikeChan)
            spike_per_trial = ComputeSpikeRates(data, num_trials, VSTIM_ON_CD, TRIAL_END_CD, 0, 0);
            h_fr = findall(gcbf,'Tag','FiringRate'); cla(h_fr);
            plot(h_fr,1:num_trials,spike_per_trial(SpikeChan,:),'.r'); hold(h_fr,'on');
            plot(h_fr,find(select_trials),spike_per_trial(SpikeChan,select_trials),'.-g');
            axis(h_fr,'tight');
            set(h_fr,'xlim',[1 num_trials]);
        end
    end
    
end

if JustEditIt == 2  % Just plot firing rates. HH20160425
    return;
end

if (~isempty(data.eye_data))
    %add the compute eye pos function here
    data.eye_positions = ComputeMeanEyePos(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
    %this next function checks to see if there is an eye calibration parameter file.  IF it exists,
    %the function will load in the parameters and compute the calibrated eye positions, and then store
    %them back in the 'data' structure for use elsewhere.  Added 12/6/00 by GCD
    data.eye_calib_done = 0;  % set this flag to zero; it can be used elsewhere to check if there are calibrated signals
    %    [caldata, doneflag] = LoadEyeCalibration(data, PATH, FILE);
    [caldata, doneflag] = LoadEyeCalibration_NonLin(data, PATH, FILE);
    data.eye_positions_calibrated = caldata;
    data.eye_calib_done = doneflag;
end

switch(Protocol)		%call a .m file that contains protocol-specific analysis routines

    %%%%%%%%%%%%%% Frequently used %%%%%%%%%%%%%%%%%%
    case AZIMUTH_TUNING_1D
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
    case DIRECTION_TUNING_3D
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
    case HEADING_DISCRIM
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
    case DELAYED_SACCADE
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
    case MEMORY_SACCADE
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, batch_flag);

        
        
         %case DIR3D_VESTIB_ACCEL
        % MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
       
    case DIRECTION_TUNING
        DirectionTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case SPEED_TUNING
        SpeedTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case HDISP_TUNING
        HDispTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case DEPTH_DISCRIM
        DepthDiscrim_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case DIRECTION_DISCRIM
        DirDiscrim_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case HDISP_GRAD_TUNING
        HGradTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SIZE_TUNING
        SizeTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURROUND_MAPPING
        SurroundMapping_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURROUND_TUNING
        SurroundTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case RELATIVE_DISPARITY
        RelativeDisparity_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case STEREOACUITY
        Stereoacuity_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case EYE_CALIBRATION
        EyeCalib_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case RF_MAPPING
        RFMapping_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case TRANS_RELATIVE_DISPARITY
        TransRelativeDisparity_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SIM_DIST_DISP_ONLY
        SimDistDispOnly_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case DEPTH_DISCRIM_NOVAR
        DepthDiscrim_NOVAR_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SIM_DIST_VERG_ONLY
        SimDistVergOnly_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SIM_DIST_DISP_VERG
        SimDistDispVerg_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SIM_DIST_CURVATURE_DISCRIM
        SimDistCurvature_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case AXIS_CUED_DIRECTION_TUNING
        AxisCuedDirecTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case CUED_DIRECTION_DISCRIM
        CuedDirectionDiscrim_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case ORIENT_CUE_DIRECTION_TUNING
        OrientCueDirectionTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
        % case DELAYED_SACCADE
        %     DelayedSaccade_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case MOVING_TARGET
        MovingTarget_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case VAR_REWARD_DIRECTION_DISCRIM
        VarRewardDirectionDiscrim_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    
    case HEADING_DISCRIM_LIP
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE); % added by Tunde 06/20/08
    case ADAPT_HEADING_DISCRIM
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE); % added by Tunde 11/06/09
    case HEADING_DISCRIM_FIXONLY
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case DIR3D_VARY_FIXATION
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZ_1D_VARY_FIXATION
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case EL_1D_VARY_FIXATION
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZ_1D_VARY_FIX_HEADEYE
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case DIR2D_CUE_CONFLICT
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case DIRECTION_REVCORR
        RevCorr_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case ROTATION_TUNING_3D
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case RVOR_PURSUIT
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case TILT_TRANSLATION
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE); %added by AHC
    case HEADING_DISCRIM_2I
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case SINUSOID
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case PURSUIT_HEADING_DISCRIM
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZIMUTH_TUNING_1D_TRAP
        MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZ_TUNUNG_W_PURSUIT_VISUAL
        Pursuit_Visual_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZ_TUNUNG_W_PURSUIT_COMBINE
        Pursuit_Combine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case AZ_TUNUNG_W_PURSUIT
        Pursuit_VisualCombine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
        %-------------------------------------------
        %Analyses for Surface Orientation Stuff, JDN 01/22/04
        %-------------------------------------------
    case SURF_DEPTH_TUNING
        SurfDepth_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURF_SPEED_TUNING
        SurfSpd_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURF_RF_MAPPING
        SurfRFMap_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURF_DIRECTION_TUNING
        SurfDirectTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURF_SIZE_TUNING
        SurfSizeTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SURF_TUNING
        SurfTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
    case SLANT_TUNING
        SurfTuning_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);
        
        %-------------------------------------------
        %Analyses for Motion parallax protocols - JWN 09/16/04
        %-------------------------------------------
    case MOTION_PARALLAX_FIX
        MPFix_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
        
        %insert analyses for binding protocols
    case FIXATION
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case FIX_1_23
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case FIX_1_23_45
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case FIX_VARY_HISTORY
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case FIX_VARY_BACKCOLOR
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
        
    case BIND_DIR_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_SPATFREQ_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_TEMPFREQ_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_BAR_DIR_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_BAR_SPEED_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_BAR_HDISP_TUNING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_RF_MAPPING
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin,StopEventBin, PATH, FILE);
    case BIND_DIR_DISC
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_BACKCOLOR
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_BACKHDISP
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_OBJ_LUM
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_OBJ_WIDTH
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_AP_RLUM
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_OBJ_POSITION
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_DOT_DENSITY
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_COHER_DISC
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case FIX_TRANS_QUAD
        CircTrans_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    case BIND_CLOSURE_DISC
        Binding_Analyses(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
        
    otherwise
        set(ListHandle, 'String','The Protocol # is not known.  Check the list of Protocol codes in TEMPO_Defs.m');
end

return;