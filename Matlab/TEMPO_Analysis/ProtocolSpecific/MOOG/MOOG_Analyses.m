function MOOG_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, batch_flag);

if nargin <15  %  no batch_flag. HH20140510
    batch_flag = [];
end

global JustEditIt;

switch(Analysis{1})
    
    %%%%%%%%%%%%%% Frequently used by HH%%%%%%%%%%%%%%%%%%
    
    % Tuning related
    case 'Plot Tuning Azimuth_HH'
        if JustEditIt == 1
            edit DirectionTuningPlot_1D_HH;
        else
            DirectionTuningPlot_1D_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'DirectionTuningPlot_1D_pairwiseunits_HH'
        if JustEditIt == 1
            edit DirectionTuningPlot_1D_pairwiseunits_HH;
        else
            DirectionTuningPlot_1D_pairwiseunits_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        end
    case 'Plot Tuning Azimuth PSTH_HH',
        if JustEditIt == 1
            edit Azimuth_PSTH_HH;
        else
            Azimuth_PSTH_HH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
        end
    case 'Plot Tuning Azimuth PSTH (for ilPPC)_HH'
        Azimuth_PSTH_HH_for_ilPPC(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);        
    case 'Plot Tuning Azimuth PSTH (for ilPPC, 2 coh)_HH'
        Azimuth_PSTH_HH_for_ilPPC_2coh(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);        
    case 'Plot Tuning Surface_HH'
        if JustEditIt == 1
            edit DirectionTuningPlot_3D_HH;
        else
            DirectionTuningPlot_3D_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'DirectionTuningPlot_3D_pairwiseunits_HH'
        if JustEditIt == 1
            edit DirectionTuningPlot_3D_pairwiseunits_HH;
        else
            DirectionTuningPlot_3D_pairwiseunits_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        end
            
    % Heading discrimination related
    case 'Plot Psychometric_HH'
        if JustEditIt == 1
            edit Psychometric_HH;
        else
            Psychometric_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'Plot Psychometric_HH_dt_yuchen'
        if JustEditIt == 1
            edit Psychometric_HH_dt_yuchen;
        else
            Psychometric_HH_dt_yuchen(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'Plot CP_HH'
        if JustEditIt == 1
            edit HeadingDis_cum_HH;
        else
            HeadingDis_cum_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'Plot CP Distribution_HH'
        if JustEditIt == 1
            edit Heading_CP_Distrib_HH;
        else
            Heading_CP_Distrib_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        end
    case 'Plot CP_shiftwindow_HH'
        if JustEditIt == 1
            edit HeadingDis_cum_shiftwindow_HH;
        else
            HeadingDis_cum_shiftwindow_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'Plot HeadingDiscimination_PSTH_HH'
        if JustEditIt == 1
            edit HeadingDis_cum_PSTH_HH;
        else
            HeadingDis_cum_PSTH_HH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);
        end
    case 'Plot HeadingDiscimination_PSTH_HH_dt_yuchen'
        if JustEditIt == 1
            edit HeadingDis_cum_PSTH_HH_dt_yuchen;
        else
            HeadingDis_cum_PSTH_HH_dt_yuchen(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);
        end        
    case 'Plot HeadingDiscimination_PSTH_HH_dt_yuchen_patch_for_HD'
        if JustEditIt == 1
            edit HeadingDis_cum_PSTH_HH_dt_yuchen_patch_for_HD;
        else
            HeadingDis_cum_PSTH_HH_dt_yuchen_patch_for_HD(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);
        end        
    case 'HeadingDis_cum_pairwiseunits_HH'
        if JustEditIt == 1
            edit HeadingDis_cum_pairwiseunits_HH;
        else
            HeadingDis_cum_pairwiseunits_HH(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        end
    case 'Plot Microstim_HH'
        if JustEditIt == 1
            edit Heading_microstim_HH;
        else
            Heading_microstim_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
        end
        
        % Others
    case 'Delayed Saccade Analysis_HH'
        if JustEditIt == 1
            edit Delayed_Saccade_Analysis_HH;
        else
            Delayed_Saccade_Analysis_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end
    case 'Memory Saccade Analysis_HH'
        if JustEditIt == 1
            edit Memory_Saccade_Analysis_HH;
        else
            Memory_Saccade_Analysis_HH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
        end

        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        
    case 'Plot Tuning Azimuth'
        DirectionTuningPlot_1D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot Tuning Surface_yong'
        DirectionTuningPlot_3D_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot Tuning Surface'
        DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot Psychometric'
        Psychometric(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot Psycho_neuro_cum'
        HeadingDis_cum(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot Microstim'
        Heading_microstim(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, batch_flag);
    case 'Delayed Saccade Analysis (Yong)'
        Delayed_Saccade_Analysis_Yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'shiftwindow'
        HeadingDis_cum_shiftwindow(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
    case 'Plot HeadingDiscimination_PSTH'
        HeadingDis_cum_PSTH(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);

        
    case 'Plot Tuning Surface (Aihua)',
        DirectionTuningPlot_3D_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Mean Variance Firing Output (Tunde)',
        Mean_Var_Firing_Output(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface (Tanya)',
        DirectionTuningPlot_3D_Tanya(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning SurfaceLFP_yong'
        DirectionTuningPlot_3DLFP_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface (Sheng)'
        DirectionTuningPlot_3D_sheng(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_LFP_xiaodong'
        DirectionTuningPlot_3D_LFP_xiaodong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Plot Lambert Tuning Surface'
        Lambert_DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'TuningPlot'
        TuningPlot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin,PATH, FILE);
    case 'Plot Tuning Azimuth_Adhira'
        DirectionTuningPlot_1D_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);   
    case 'Direction/Disparity Tuning (Aihua)'
        Direction_Disparity_Tuning_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);         
    case 'Plot Tuning Azimuth PSTH',
        Azimuth_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'DirectionTuningPlot_3D_pairwiseunits_yong'
        DirectionTuningPlot_3D_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_3D_pairwiseunits_sam'
        DirectionTuningPlot_3D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_eyetrace'
        DirectionTuningPlot_1D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_firingrate'
        DirectionTuningPlot_1Dfiringrate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_firingrate (Tunde)'
        DirectionTuningPlot_1Dfiringrate_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'DirectionTuningPlot_1D_pairwiseunits_yong'
        DirectionTuningPlot_1D_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_pairwiseunits_sam'
        DirectionTuningPlot_1D_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_ves'
        Directiontuningplot_ves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Event Times'
        EventTimes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot PSTH'
        MOOG_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH (Tanya)'
        MOOG_PSTH_Tanya(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_yong'
        PSTH_Translation_yong(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'MOOG_PSTH_xiongjie'
        MOOG_PSTH_xiongjie(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk'
        PSTH_Anuk(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis'
        PSTH_Corr_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis_fix'
        PSTH_Corr_Analysis_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk_fix'
        PSTH_Anuk_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run Frequency Analysis'
        Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run Frequency Analysis alternate'
        Frequency_Analysis_alternate(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot PSTH_Anuk_HTI_fix'
        PSTH_Anuk_HTI_fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Spike Rasters'
        PlotRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, Protocol);
    case 'Plot Psycho_neuro_cum_Adhira'
        HeadingDis_cum_Adhira(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cum (Sheng)'
        HeadingDis_cum_Sheng(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Heading_tuning'
        HeadingDis_fixation(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Reaction Time Task)'
        Psychometric_RT(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Reaction Time Analysis (LIP_Adhira)'
        Reaction_Time_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Psychometric Conflict'
        Psychometric_1I_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric (Adaptation)'
        Psychometric_Adaptation(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);                
    case 'Plot Psychometric (Adam)'
        Psychometric_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Psychometric (Adaptation, Adam)'
        Psychometric_Adaptation_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Psychometric_pursuit'
        Psychometric_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot PURSUIT_HEADING_eyemovement'
        PURSUIT_HEADING_eyemovement(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Accelerometer_cum'
        Accelerometer_cum(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Accelerometer_cum (Tunde)'
        Accelerometer_cum(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot CP Distribution'
        Heading_CP_Distrib(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Calculate HeadingDiscimination_DFT'
        HeadingDis_cumDFT(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Calculate HeadingDiscimination_Corr'
        HeadingDis_cumCorr(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'HeadingDis_cum_eyetrace'
        HeadingDis_cum_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_eyetrace_staircase'
        HeadingDis_cum_eyetrace_staircase(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cumVarience'
        HeadingDis_cumVarience(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psycho_neuro_cumTimecourse'
        HeadingDis_cumTimecourse(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot session_yong'
        HeadingDis_cum_session(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_pairwiseunits_yong'
        HeadingDis_cum_pairwiseunits_yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'HeadingDis_cum_pairwiseunits_sam'
        HeadingDis_cum_pairwiseunits_sam(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_Fix_sam'
        DirectionTuningPlot_Fix_SUSU(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'DirectionTuningPlot_1D_Az_VF_sam'
        DirectionTuningPlot_1D_Az_VF_SUSU(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'LessSampleFit'
        LessSampleFit(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Eye trace'
        DirectionTuning3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Eye trace (tunde)'
        DirectionTuning3D_eyetrace_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);            
    case 'Plot Fixation Tuning Surface'
        DirectionTuningPlot_Fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Fixation Tuning Surface (Tunde)'
        DirectionTuningPlot_Fix_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Plot Fixation Tuning Surface (Aihua)'
        DirectionTuningPlot_Fix_aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);     
    case 'Wrapped Gaussian Fit (Tunde)'
        WrappedGaussianFit_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Lambert Fixation Tuning Surface',
        Lambert_DirectionTuningPlot_3D_fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Fixation PSTH'
        MOOG_PSTH_Fix(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Output Data for Curve Fitting'
        DirectionTuningPlot_Curvefit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_conflict'
        Direction2d_cue_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning Surface_conflict (Aihua)'
        Direction2d_cue_conflict_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Fit_cue_conflict_2D'
        Direction2d_cue_conflict_fit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Info_cue_conflict_2D'
        Direction2d_cue_conflict_info(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Basic_cue_conflict_2D'
        Direction2d_cue_conflict_basic(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tuning_yong'
        Direction2d_cue_conflict_basic_yong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Fit Optic Flow Tuning (Zack)'
        Heading_CurveFit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Fit Optic Flow Tuning_Rot (Zack)' % for rotation
        Heading_CurveFit_Rot(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Azimuth'
        DirectionTuningPlot_1D_Az_VF(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Azimuth (Tunde)'
        DirectionTuningPlot_1D_Az_VF_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation Tuning Curves_Elevation'
        DirectionTuningPlot_1D_El_VF(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot 1D Fixation PSTH'
        MOOG_PSTH_Fix_1D(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot 1D Fixation Tuning Curves_Azimuth Head-Eye'
        DirectionTuningPlot_1D_Az_VF_headeye(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Rotation Tuning 3D'
        Rotation3Dtuning(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Plot Rotation Tuning 3D (Aihua)'
        Rotation3Dtuning_Aihua(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);       
    case 'Plot Rotation Tuning 3D (Adhira)'
        Rotation3Dtuning_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);  
    case 'Plot Rotation Tuning 3D (Sheng)'
        Rotation3Dtuning_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);                   
    case 'Plot Lambert Rotation Tuning 3D'
        Lambert_Rotation3Dtuning(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case '2D Frequency_Analysis'
        Frequency_Analysis_2d(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Rotation Frequency Analysis'
        Rotation_Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH'
        Rotation_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH (Tanya)'
        Rotation_PSTH_Tanya(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);    
     case 'Plot Rotation PSTH_xiongjie'
        Rotation_PSTH_xiongjie(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation PSTH_yong'
        PSTH_yong(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Rotation Eye Trace (Katsu)'
        Rotation3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Rotation Eye Trace (Tunde)'
        Rotation3D_eyetrace_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Firing rate'
        FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Firing rate (Tunde)'
        FiringRate_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Output ROT_Firing rate'
        Rotation_FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot  TILT_TRANSLATION PSTH'
        TiltTrans_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);  %added by AHC
    case 'Plot Tilt/Trans Direction tuning'
        DirectionTuning2D_TiltTrans(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Pursuit Direction Tuning'
        DirectionTuning2D_pursuit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Pursuit PSTH'
        DirectionTuning2D_pursuit_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);%added by KT
    case 'MU activity'
        %DirectionTuningPlot_3D_munit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        MunitAnalysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'MU activity (Xiaodong)'
        %DirectionTuningPlot_3D_munit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
        MunitAnalysis_Xiaodong(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Plot PSTH_Anuk_rotation'
        PSTH_Anuk_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Corr_Analysis_rotation'
        PSTH_Corr_Analysis_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Run PSTH_Gaussfit'
        %PSTH_Anuk_rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
        PSTH_Gaussfit(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Psychometric 2I'
        Psychometric_2I(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Psychometric 2I Conflict'
        Psychometric_2I_conflict(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Sinusoid Analysis'
        Sinusoid_Analysis(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Accelerlometer study'
        Accelerometer_study(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Tilt/Trans Eye Trace'
        TiltTrans_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Horizontal only (Katsu)'
        Horizontal_only(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Plot Pursuit Eye Trace'
        Pursuit_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Rotation Horizontal only (Katsu)'
        Rotation_Horizontal_only(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Rotation_Frontal (Katsu)'
        Rotation_Frontal(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);
    case 'Save 4direct Pursuit Eye Trace'
        Pursuit_eachCell_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot  First1sec_PSTH'
        First1sec_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);  %added by Katsu
    case 'FanoFactor (Katsu)'
        FanoFactor(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Rotation3D_DDIperm (Katsu)'
        Rotation3D_DDIperm(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Direction3D_DDIperm (Katsu)'
        Direction3D_DDIperm(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Temporal_Analysis_Translation (Katsu)'
        Temporal_Analysis_Translation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Temporal_Analysis_Rotation (Katsu)'
        Temporal_Analysis_Rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Compute Vertical Slice'
        VerticalSlice(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Lambert Rotation Tuning 3D (Katsu)'
        Katsu_Lambert_plotter(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Output Rotation PSTH (Katsu)'
        Rotation_PSTH_output(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Analysis);
    case 'Output Translation PSTH (Katsu)'
        Translation_PSTH_output(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Analysis);
    case 'Frontal parralel plane'
        Frontal_parralel_plane(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Temporal Analysis (cah)'
        HeadingDis_temporal_cah(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot_PSTH_Michael'
        MOOG_PSTH_Michael(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'DirectionTuningPlot_1D (Tunde)'
        DirectionTuningPlot_1D_Az_VF_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuning_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_Tuning_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuningcoherence_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_Coherence_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case  'Plot 1Dtuningtimecourse_trapezoid (Yong)'
        AZIMUTH_TUNING_1D_TRAP_timecourse_Yong(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Plot Accelerometer_cum_1D (Jing)'
        Accelerometer_cum_1D_Jing(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Delayed Saccade Analysis (Adhira)'
        Delayed_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);               
    case 'Delayed Saccade Analysis (Adam)'
        Delayed_Saccade_Analysis_Adam(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
    case 'Memory Saccade Analysis (Adhira)'
        Memory_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   case 'Memory Saccade Analysis (Gu)'
        Memory_Saccade_Analysis_Gu(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

    case 'Galvo Pursuit'
        Delayed_Saccade_Analysis_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    case 'Pursuit Visual (Adhira)'
        Pursuit_Visual_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);   
    case 'Pursuit-Combine (Adhira)'
        Pursuit_Combine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
    case 'Pursuit -Visual+Combine (Adhira)'
        Pursuit_VisualCombine_Adhira(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);        
end
    
return;