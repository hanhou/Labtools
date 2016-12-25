
% pwatkins - START fitting specific options

fitting_stims = [CURVEFIT_OCCULAR_STIM];  % which stim types to fit to


fitting_models = [@Curvefit_head_centered_7p_asp_ng ];
fitting_bounds = [@Curvefit_head_centered_bounds_7p_asp_ng ];
fitting_tuning_models = [@Curvefit_cos_tuning_7p_asp ];
%fitting_models = [@Curvefit_head_centered_7p_asp_ng @Curvefit_head_centered_8p_asp_ng];
%fitting_bounds = [@Curvefit_head_centered_bounds_7p_asp_ng @Curvefit_head_centered_bounds_8p_asp_ng];
%fitting_tuning_models = [@Curvefit_cos_tuning_7p_nof_asp @Curvefit_cos_tuning_8p_asp ];
%fitting_models = [@Curvefit_head_centered_ng];
%fitting_bounds = [@Curvefit_head_centered_bounds_ng ];
%fitting_tuning_models = [@Curvefit_cos_tuning_7p ];



minimization_iterations = 100;  % number of times to iterate a complete minimization for each model fit
fitting_options = optimset('MaxFunEvals',10000,'MaxIter',1000,...  % minimization options for fmincon
    'LargeScale','off','MaxSQPIter',5000);

% which plots to generate
plot_contour = 1;
plot_R2_export = 1;  % export only, plot somewhere else
plot_param_export = 1;  % export only, plot somewhere else
plot_R2_distrib_export = 0;  % export only, plot somewhere else
plot_R2_distrib = 0;

% options for contour plots
plot_clabels = 0;   % whether to plot data labels in contour plots
plot_str_len = 70;  % maximum string length for special statistic text
start_figure_contour = 2;   % figure number to start at
plot_contour_num_compares = 1;  % number of fits to compare with next fit
plot_contour_compares_step = 2;  % number of fits to compare with next fit
save_contour_figures = 0;
print_contour_figures = 0;
printer_contour_figures = '-PTektronix Phaser 750N';

% options for bootstrap r^2 distrib plots
start_figure_R2_distrib = 50;
plot_R2_distrib_num_bins = 10;
save_R2_distrib_figures = 0;
print_R2_distrib_figures = 0;
printer_R2_distrib_figures = '-PHP LaserJet 4200dtnsl';

% options for bootstrap stats
compute_stat_bootstrap = 0;
bootstrap_num_samples = 10000;

% option to only fit single gaze angle
fit_single_gaze_angle = 0;

use_backdoor_load = 1;
use_backdoor_fit = 1;
build_backdoor_load = 0;
build_backdoor_fit = 0;

% which phases of the fitting to run, dependent in sequence.
run_load = 1;
run_fit = 1;
run_compute_stat = 1;
run_plot = 1;

% backdoor preference?
if use_backdoor_load | use_backdoor_fit | build_backdoor_load | build_backdoor_fit | ...
        save_contour_figures | save_R2_distrib_figures
    backdoor_dir = fullfile('work','tempo_backdoor');
    % create the backdoor directory, incase it does not exist.
    % for some stupid reason matlab can not create an absolute directory.
    tmp = pwd;
    cd(matlabroot);
    [s,mess,messid] = mkdir(backdoor_dir);
    cd(tmp);
    backdoor_dir = fullfile(matlabroot, backdoor_dir);
    backdoor_load_ext = '_load.mat';
    backdoor_fit_file = fullfile(backdoor_dir, 'results.txt');
    save_figure_ext = '.fig';
    save_figure_format = 'fig';
end

% END fitting specific options
