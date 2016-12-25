function result = CP_HH(headings,choices,spike_counts,CP_p_method, IF_FIGURE)
% result = CP_HH(headings,choices,spike_counts, {CP_p_method, {IF_FIGURE}})
%   Calculate choice probability and neuro/psycho- metric functions.
%   This code is more concise, clear, and reusable than the previous
%   messy versions which I cannot bear anymore...
%
%   CP_p_method = 0, no p value returned
%               > 0, permutation method with CP_p_method times, and also return the t-test p value for comparison
%               = -1, use two sample t-test p value (much faster)
%
%   @HH20140925

% tic
% Parallel computing
if matlabpool('size') == 0
    try 
        matlabpool; 
    catch
    end
end  

% if  isempty(gcp('nocreate'))
%     try parpool(4); catch ; end
% end  
    
if nargin < 4
    CP_p_method = 1000; % Default: permutate 1000 times
    IF_FIGURE = 0;
elseif nargin < 5
    IF_FIGURE = 0;
end

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%
method = 0; % 0: Maximum likelihood; 1: Square error
tolerance = 10; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headings = headings(:);
choices = choices(:);
spike_counts = spike_counts(:);

LEFT = 1;
RIGHT = 2;

unique_heading = munique(headings);

% Sometimes I just reuse the codes in CP_HH to calculate auROC and p_value
% for other purposes (see ChoicePreference and ModalityPreference in
% HeadingDis_cum_PSTH_HH), where the input heading, choice, etc., are therefore called "fake". @hh20150418
if_fake = ~(length(unique_heading) > 1 && abs(sum(unique_heading)) < 0.0001);  % More than one heading and symmetric

%% Psychometric function (Only when we need a figure)

if IF_FIGURE
    
    % Psycho function
    for hh = 1:length(unique_heading)
        num_headings(hh,1) = sum(headings == unique_heading(hh));
        rightward_prop(hh,1) = sum(choices(headings == unique_heading(hh)) == RIGHT)/ sum(headings == unique_heading(hh));
    end
    
    % Fitting
    [Psy_bias,Psy_thres] = cum_gaussfit_max1([unique_heading, rightward_prop, num_headings], method, 0);
    [Psy_bias_tol,Psy_thres_tol] = cum_gaussfit_max1([unique_heading, rightward_prop, num_headings], method, tolerance);
   
    % Saving
    result.Psy_func = [unique_heading, rightward_prop, num_headings];
    result.Psy_para = [Psy_bias,Psy_thres];
    result.Psy_para_tol = [Psy_bias_tol,Psy_thres_tol];
    
    figure(222);
    plot(unique_heading,rightward_prop,'ko'); hold on
    xx = min(unique_heading):0.01:max(unique_heading);
    h1 = plot(xx,cum_gaussfit(result.Psy_para,xx),'k');
end

%% Neurometric function

if ~if_fake
    
    correct_trials = (sign(choices - 1.5) == sign(headings)) | (sign(headings)==0);
    
    % ---- Neuro local tuning curve ----
    for hh = 1:length(unique_heading)
        curr_heading = headings == unique_heading(hh);
        resp_mean(hh,1) = mean(spike_counts(curr_heading));
        resp_se(hh,1) = std(spike_counts(curr_heading)) / sqrt(sum(curr_heading));
        
        % Another version which only includes correct trials. HH20150403
        curr_heading_correctonly = (headings == unique_heading(hh)) & correct_trials;
        resp_mean_correctonly(hh,1) = mean(spike_counts(curr_heading_correctonly));
        resp_se_correctonly(hh,1) = std(spike_counts(curr_heading_correctonly)) / sqrt(sum(curr_heading_correctonly));
        
    end
    
    [rr,pp] = corrcoef(unique_heading, resp_mean);
    result.pref = (rr(1,2) <= 0) * LEFT + (rr(1,2) > 0) * RIGHT;
    
    result.Neu_tuning = [unique_heading, resp_mean, resp_se];
    
    result.Neu_tuning_correctonly = [unique_heading, resp_mean_correctonly, resp_se_correctonly];
    
    % Plotting
    if IF_FIGURE
        figure(223);
        errorbar(unique_heading, resp_mean, resp_se, 'o-');
        title('Local Tuning');
    end
    
    % ------ Neurometric function -------
    % 1. Compared with 0 heading
    if sum(unique_heading == 0) == 0 % If we don't have 0 heading, skip it
        rightward_prop_neuro_with0 = ones(sum(unique_heading ~= 0),1) * NaN;
    else
        for hh = 1:length(unique_heading)
            % Assuming RIGHT preferred direction first, then rocN(response to theta, response to 0)
            % should be >0.5 if theta > 0, and <0.5 if theta <0, a resonable result.
            rightward_prop_neuro_with0(hh,1) = rocN(spike_counts(headings == unique_heading(hh)),spike_counts(headings == 0));
        end
        
        % Flip neurometrtic function if the preferred direction is LEFT
        if result.pref == LEFT
            rightward_prop_neuro_with0 = 1 - rightward_prop_neuro_with0;
        end
        
        % Subtract 0 heading from neurometric function
        rightward_prop_neuro_with0 (unique_heading == 0) = [];
    end
    
    %----------------
    
    % 2. Compared with anit-neuron
    for hh = 1: sum(unique_heading ~=0) / 2 % Neurometric function with anti-neuron is symmetric, so we just calculate the first half.
        % Assuming RIGHT preferred direction first, then rocN(response to theta, response to -theta)
        % should be >0.5 if theta > 0, and <0.5 if theta <0, a resonable result.
        rightward_prop_neuro_anti(hh,1) = rocN(spike_counts(headings == unique_heading(hh)),spike_counts(headings == - unique_heading(hh)));
    end
    
    % Flip if the preferred direction is LEFT
    if result.pref == LEFT
        rightward_prop_neuro_anti = 1 - rightward_prop_neuro_anti;
    end
    
    % Mirror to the other half
    rightward_prop_neuro_anti = [rightward_prop_neuro_anti; 1 - flipud(rightward_prop_neuro_anti)];
    
    
    %----------------
    % Fitting
    %  [Neu_bias_with0,Neu_thres_with0] = cum_gaussfit_max1([unique_heading(unique_heading~=0), rightward_prop_neuro_with0]);
    [Neu_bias_with0,Neu_thres_with0] = deal(NaN,NaN);
    [Neu_bias_anti,Neu_thres_anti] = cum_gaussfit_max1([unique_heading(unique_heading~=0), rightward_prop_neuro_anti],method,0); % No tolerance for neurometric curve
    
    % Negative and positive infinite value means flat tuning
    if Neu_thres_with0 < 0 || Neu_thres_with0> 300
        Neu_thres_with0 = 300;
    end
    
    if Neu_thres_anti < 0 || Neu_thres_anti > 300
        Neu_thres_anti = 300;
    end
    
    % Saving
    result.Neu_func_with0 = [unique_heading(unique_heading~=0), rightward_prop_neuro_with0];
    result.Neu_para_with0 = [Neu_bias_with0,Neu_thres_with0];
    
    result.Neu_func_anti = [unique_heading(unique_heading~=0), rightward_prop_neuro_anti];
    result.Neu_para_anti = [Neu_bias_anti,Neu_thres_anti];
    
    
    
    % Plotting
    if IF_FIGURE
        figure(222);  hold on;
        xx = min(unique_heading):0.01:max(unique_heading);
        h2 = plot(result.Neu_func_with0(:,1),result.Neu_func_with0(:,2),'rs');
        plot(xx,cum_gaussfit(result.Neu_para_with0,xx),'r-');
        
        h3 = plot(result.Neu_func_anti(:,1),result.Neu_func_anti(:,2),'m^');
        plot(xx,cum_gaussfit(result.Neu_para_anti,xx),'m-');
        
        legend([h1 h2 h3],{'Psycho','Neuro with 0','Neuro anti'},'Location','Best');
    end
     
else % If faked, just let the preferred to be RIGHT
    result.pref = RIGHT;
end
%% Choice probability

% ------  Group data according to choices and calucate CP for all headings ------
spike_counts_allheadings_grouped = cell(length(unique_heading),2);

for hh = 1:length(unique_heading)
    curr_heading = headings == unique_heading(hh);
    
    % We only consider headings that have at least 3 LEFT & RIGHT choices.
    if (sum(curr_heading & (choices == LEFT)) >= 3) && (sum(curr_heading & (choices == RIGHT)) >= 3)
        % Raw spike counts for each headings
        spike_counts_allheadings_grouped{hh,LEFT} = spike_counts(curr_heading & (choices == LEFT));
        spike_counts_allheadings_grouped{hh,RIGHT} = spike_counts(curr_heading & (choices == RIGHT));
        
        % Calculate CP for current heading: rocN(pref,null)
        result.CP_allheadings(hh) = rocN(spike_counts_allheadings_grouped{hh, result.pref}, ...
            spike_counts_allheadings_grouped{hh, LEFT + RIGHT - result.pref});
    else
        result.CP_allheadings(hh) = NaN;
    end
    
end

result.spike_counts_allheadings_grouped = spike_counts_allheadings_grouped;

% Automatically change the "dead-ahead" if the bias is huge
% if abs(Psy_bias) > 0.5
%     [~,h_0]=min(abs(Psy_bias-unique_heading)); % HH20130905
%     %      [~,h_0]=min(abs(rightward_prop-0.5)); % HH20140510
% else
%     h_0 = find(unique_heading >= 0,1);
% end
h_0 = find(unique_heading >= 0,1);

% CP for "0 heading"
result.dead_ahead = unique_heading(h_0);
result.CP_0 = result.CP_allheadings(h_0);

% ----- Calculate Grand CP -----

if length(unique_heading)>1
    
    % Z-score
    spike_zscored = ones(size(spike_counts)) * NaN;
    spike_zscored_grand_grouped = cell(1,2);
    
    for hh = 1:length(unique_heading)
        curr_heading = headings == unique_heading(hh);
        
        % We only consider headings that have at least 3 LEFT & RIGHT choices.
        if (sum(curr_heading & (choices == LEFT)) >= 3) && (sum(curr_heading & (choices == RIGHT)) >= 3)
            % Z-score current heading
            z_dist = spike_counts(curr_heading);
            z_dist = (z_dist - mean(z_dist))/std(z_dist);
            spike_zscored(curr_heading) = z_dist;
        end
    end
    
    % Group z-scored firing rates for all headings and calculate grand CP
    spike_zscored_grand_grouped{LEFT} = spike_zscored(choices == LEFT & ~isnan(spike_zscored));
    spike_zscored_grand_grouped{RIGHT} = spike_zscored(choices == RIGHT & ~isnan(spike_zscored));
    result.CP_grand = rocN(spike_zscored_grand_grouped{result.pref}, spike_zscored_grand_grouped{LEFT + RIGHT - result.pref});
    
    result.spike_zscored = spike_zscored;
    result.spike_zscored_grand_grouped = spike_zscored_grand_grouped;
    
end

% ---- P value for CPs -----

% Method 1: Permutation test

if CP_p_method > 0

    % P value for CP at "0 heading"
    if ~isnan(result.CP_0)
        spike_counts_0_all_choices = cat(1,spike_counts_allheadings_grouped{h_0,:});             % Data to be permuted
        CP_0_perm = nan(CP_p_method,1);        % Pre-allocation
        n_pref = length(spike_counts_allheadings_grouped{h_0,result.pref});
        
        parfor i = 1 : CP_p_method   % For each permutation (total times: CP_p_method)
            % Permuted CP_0
            rand_choice_pref = randperm(length(spike_counts_0_all_choices),n_pref);
            rand_choice_null = setdiff(1:length(spike_counts_0_all_choices),rand_choice_pref);
            CP_0_perm(i) = rocN(spike_counts_0_all_choices(rand_choice_pref),spike_counts_0_all_choices(rand_choice_null));
        end
        
%         data_GPU = gpuArray(spike_counts_0_all_choices);
%         n_pref = length(spike_counts_allheadings_grouped{h_0,result.pref});
%         CP_0_perm = arrayfun(@PermutedCP, num2cell(repmat(data_GPU,1,1),1),repmat(n_pref,1,1),'UniformOutput',0);
        
        % Two-tail test
        if result.CP_0 > 0.5
            result.CP_0_p_perm = 2 * (1 - sum(result.CP_0 > CP_0_perm) / length(CP_0_perm));
        else
            result.CP_0_p_perm = 2 * sum(result.CP_0 > CP_0_perm) / length(CP_0_perm);
        end
        
    else
        result.CP_0_p_perm = NaN;
    end
    
    if length(unique_heading)>1
        
        % P value for grand CP (same as above)
        if ~isnan(result.CP_grand)
            spike_zscored_all_choices = cat(1,spike_zscored_grand_grouped{:});  % Data to be permuted
            CP_grand_perm = zeros(CP_p_method,1);  % Pre-allocation
            
            parfor i = 1 : CP_p_method   % For each permutation (total times: CP_p_method)
                % Permuted Grand CP
                rand_choice_pref = randperm(length(spike_zscored_all_choices),length(spike_zscored_grand_grouped{result.pref}));
                rand_choice_null = setdiff(1:length(spike_zscored_all_choices),rand_choice_pref);
                CP_grand_perm(i) = rocN(spike_zscored_all_choices(rand_choice_pref),spike_zscored_all_choices(rand_choice_null));
            end
            
            % Two-tail test
            if result.CP_grand > 0.5
                result.CP_grand_p_perm = 2 * (1 - sum(result.CP_grand > CP_grand_perm) / length(CP_grand_perm));
            else
                result.CP_grand_p_perm = 2 * sum(result.CP_grand > CP_grand_perm) / length(CP_grand_perm);
            end
        else
            result.CP_grand_p_perm = NaN;
        end
        
    end
else % CP_p_method <=0
        result.CP_grand_p_perm = NaN;
end

% Method 2: Two sample t-test (parametric, but much faster)
if CP_p_method ~= 0
    [~, result.CP_0_p_ttest] = ttest2(spike_counts_allheadings_grouped{h_0,LEFT},spike_counts_allheadings_grouped{h_0,RIGHT});
    
    if length(unique_heading)>1
        [~, result.CP_grand_p_ttest] = ttest2(spike_zscored_grand_grouped{LEFT},spike_zscored_grand_grouped{RIGHT});
    end
end

% toc


function CP_perm = PermutedCP(spike_counts_0_all_choices,n_pref)

rates = spike_counts_0_all_choices{:};

rand_choice_pref = randperm(length(rates),n_pref);
rand_choice_null = setdiff(1:length(rates),rand_choice_pref);
CP_perm = rocN_GPU(rates(rand_choice_pref),rates(rand_choice_null));
