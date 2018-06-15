function function_handles = Group_Psychometric(XlsData,if_tolerance)

%% Get data
monkey = get(findall(gcbf,'tag','monkey'),'value');

num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;

colors = {'b','r','g'};

% ==== Get data ====

if monkey == 2 % Messi
    mat_address = 'Z:\Data\Tempo\Batch\20150725_Psycho_m10\';
%     mat_address = 'Z:\Data\Tempo\Batch\20150105_MessiTraining\';
    mask_all = (strcmp(txt(:,header.Protocol),'Training')| strcmp(txt(:,header.Protocol),'HD') | strcmp(txt(:,header.Protocol),'HD_dt')) & (num(:,header.Monkey) == 10) & ~strcmp(txt(:,header.Note),'nnn');
elseif monkey == 1 % Polo
    mat_address = 'Z:\Data\Tempo\Batch\20150725_Psycho_m5\';
%     mat_address = 'Z:\Data\Tempo\Batch\20150125_PoloTraining\';
    mask_all = (strcmp(txt(:,header.Protocol),'Training') | strcmp(txt(:,header.Protocol),'HD')) & (num(:,header.Monkey) == 5) & ~strcmp(txt(:,header.Note),'nnn');
end

xls_num = num(mask_all,:);
xls_txt = txt(mask_all,:);
xls_raw = raw(mask_all,:);

% Exclude duplications (more than one cells for one psychometric session)
[~,unique_session]=unique(xls_txt(:,header.FileNo),'stable');

xls_num = xls_num(unique_session,:);
xls_txt = xls_txt(unique_session,:);
xls_raw = xls_raw(unique_session,:);

% cd(mat_address);

% group_result(sum(mask_all)).raw = [];

nSession = size(xls_raw,1);
real_Session = [];

% Load files
for i = 1:nSession
    % Find .mat file
    try
        fileName = dir([mat_address xls_txt{i,header.FileNo} '_*.mat']);
        result = load([mat_address fileName(1).name]);  % Only the first appearance
        result = result.result;
        
        if result.repetitionN >= 8 && length(result.unique_motion_coherence)==1 && (isnan(result.Thresh_psy{3})||result.Thresh_psy{3} > 0.2) % Only repN > 8 and only have one coherence (and exclude abnormal thresholds)
            
            % For dt. HH20161127
            if length(result.unique_stim_type) == 5  % Only include first three conditions
                result.correct_rate = result.correct_rate(1:3);
                result.Thresh_psy = result.Thresh_psy(1:3);
                result.Bias_psy = result.Bias_psy(1:3);
                result.Thresh_psy_tol = result.Thresh_psy_tol(1:3);
                result.Bias_psy_tol = result.Bias_psy_tol(1:3);
                result.psy_thresh_shift = result.psy_thresh_shift(1:3);
                result.psy_bias_shift = result.psy_bias_shift(1:3);
                result.psy_thresh_shift_tol = result.psy_thresh_shift_tol(1:3);
                result.psy_bias_shift_tol = result.psy_bias_shift_tol(1:3);
                result.fit_data_psycho_cum = result.fit_data_psycho_cum(1:3);
            end

            % Filling in group_result
            real_Session = [real_Session i];
            group_result(length(real_Session)) = struct(result);
        end
        
        
    catch
        disp(['Read data error:  ' mat_address xls_txt{i,header.FileNo}]);
        %         keyboard;
    end
    
end

nSession = length(real_Session);
xls_num = xls_num(real_Session,:);
xls_txt = xls_txt(real_Session,:);
xls_raw = xls_raw(real_Session,:);


% Some more preprocessings
for i = 1:nSession
    group_result(i).staircase = range(group_result(i).fit_data_psycho_cum{group_result(i).unique_stim_type(1)}(:,end)) > 1; % Staircase marker
    if strcmp(group_result(i).FILE,'m10c0r2')
        group_result(i).staircase = true;  % Special case
    end
end

% Pack data
staircase = [group_result.staircase]';
CR_session =  reshape(cell2mat([group_result.correct_rate]),3,[])';

if if_tolerance
    Thresh_session = reshape(cell2mat([group_result.Thresh_psy_tol]),3,[])';
    Bias_session = reshape(cell2mat([group_result.Bias_psy_tol]),3,[])';
    Thresh_shift_session = reshape([group_result.psy_thresh_shift_tol],3,[])';
else
    Thresh_session = reshape(cell2mat([group_result.Thresh_psy]),3,[])';
    Bias_session = reshape(cell2mat([group_result.Bias_psy]),3,[])';
    Thresh_shift_session = reshape([group_result.psy_thresh_shift],3,[])';
end

Coh_session = nan(nSession,1);
for i = 1:nSession
    if ~isempty([find(group_result(i).unique_stim_type == 2) find(group_result(i).unique_stim_type == 3)])
        Coh_session(i,1) = mode(group_result(i).raw(:,2));
    end
end


% Auxillary
tf_session = xls_num(:, header.HD_TargFirst);
glass_session = xls_num(:, header.Psy_glass);
platform_session = xls_num(:,header.Psy_platform);

%% Turn sessions into days (Average every day)
day = xls_num(:,1);
unique_day = unique(day);

nDay = length(unique_day);
Thresh_day = nan(size(unique_day,1),3);
Bias_day = nan(size(unique_day,1),3);
Coh_day = nan(size(unique_day,1),1);

prediction_ratio_interaction = nan(size(unique_day,1),2); % 1:interleaved; 2: comb_alone
        
for dd = 1:length(unique_day)
    
    thresh_today = Thresh_session(day == unique_day(dd),:);
    bias_today = Bias_session(day == unique_day(dd),:);
    coh_today = Coh_session(day == unique_day(dd),:);
    
    for condition = 1:3
        % Day-average
        Thresh_day(dd,condition) = mean(thresh_today(~isnan(thresh_today(:,condition)),condition));
        Bias_day(dd,condition) = mean(bias_today(~isnan(thresh_today(:,condition)),condition));
        
        % Auxiliary
        tf_day (dd) = mode(tf_session(day == unique_day(dd)));
        glass_day (dd) = mode(glass_session(day == unique_day(dd)));
        platform_day (dd) = mode(platform_session(day == unique_day(dd)));
    end
    
    % For interaction analysis
    sigma_pred_today = 1./sqrt(1./Thresh_day(dd,1).^2 + 1./Thresh_day(dd,2).^2);
    comb_alone_today = isnan(thresh_today(:,1)) & isnan(thresh_today(:,2)) & ~isnan(thresh_today(:,3));
    interleave_today = ~any(isnan(thresh_today),2);
    if ~isempty(interleave_today)
        prediction_ratio_interaction(dd,1) = mean(thresh_today(interleave_today,3))/sigma_pred_today;
    end
    if ~isempty(comb_alone_today)
        prediction_ratio_interaction(dd,2) = mean(thresh_today(comb_alone_today,3))/sigma_pred_today;
    end
    
    % Find coherence for each day
    Coh_day(dd,1) = mean(coh_today(~isnan(coh_today(:,1)),1));
end

interaction_pairs = fliplr(prediction_ratio_interaction(~any(isnan(prediction_ratio_interaction),2),:));

%% Manually add some data
if monkey == 2 % ~isempty(strfind(mat_address,'Messi'))
    Polo_old = load('Polo_old');
    Polo_old = Polo_old.Polo_old;
    
    Polo = Polo_old(1:min(max(size(Polo_old,1),size(Thresh_day,1)),end),[1:6]);

    % Insert Messi's records from 20150519 to 20150609. HH20150609
    Messi_insert_threshold = [5.22 4.87 nan; 5.3 nan nan; 6.5 nan nan; 5.2 nan nan; 6.6 nan nan; 4.6 nan nan; 7.5 nan nan;
                             12.5 nan nan; 8 nan nan; 7.16 nan nan; ];
    Messi_insert_bias = [-1.7 0.5 nan; -2 nan nan; -2.1 nan nan; 0.8 nan nan; -0.4 nan nan; 0.7 nan nan; 2.6 nan nan;
                        1.76 nan nan; -0.5 nan nan; 2.4 nan nan; ];
                    
    Thresh_day = [Thresh_day(1:97,:); Messi_insert_threshold; Thresh_day(98:end,:)];
    Bias_day = [Bias_day(1:97,:); Messi_insert_bias; Bias_day(98:end,:)];
    tf_day  = [tf_day(1:97) ones(1,length(Messi_insert_bias)) tf_day(98:end)];
    glass_day  = [ glass_day(1:97) zeros(1,length(Messi_insert_bias)) glass_day(98:end)];
    platform_day  = [ platform_day(1:97) 301*ones(1,length(Messi_insert_bias)) platform_day(98:end)];
    Coh_day = [Coh_day(1:97,:); nan(length(Messi_insert_bias),1); Coh_day(98:end,:)];
                    
    % Renew nDay
    nDay = size(Thresh_day,1);
    
elseif monkey == 1  % ~isempty(strfind(mat_address,'Polo')) 
    
    % Add Polo old records without htb files. HH20150413
    Polo_old = load('Polo_old');
    Polo_old = Polo_old.Polo_old;
    old_threshold = Polo_old(:,[1 3 5]);
    old_bias = Polo_old(:,[2 4 6]);
    
    Thresh_day = [old_threshold; Thresh_day];
    Bias_day = [old_bias; Bias_day];
    tf_day  = [zeros(1,length(old_bias)) tf_day];
    glass_day  = [ones(1,length(old_bias)) glass_day];
    platform_day  = [301*ones(1,length(old_bias)) platform_day];
    Coh_day = [nan(length(old_bias),1) ; Coh_day];
    
    % Renew nDay
    nDay = size(Thresh_day,1);

    % -- nSession --
    Thresh_session = [old_threshold; Thresh_session];
    Bias_session = [old_bias; Bias_session];
    tf_session  = [zeros(1,length(old_bias))'; tf_session];
    glass_session  = [ones(1,length(old_bias))'; glass_session];
    platform_session  = [301*ones(1,length(old_bias))'; platform_session];
    Coh_session = [nan(length(old_bias),1) ; Coh_session];
    
    % Renew nDay
    nSession = size(Thresh_session,1);

    
end

showPlatform = 0; % ~isempty(strfind(mat_address,'Polo'));
platforms = [301 109 103 102];
platformColor = {'k','r','g','g'};

%% ====================================== Function Handles =============================================%

function_handles = {
    'Training days', {
    'Threshold and bias', @f1p1;
    'Prediction ratio', @f1p2;
    'Average threshold', @f1p3;
    };
    
    'Training sessions', {
    'Threshold and bias', @f2p1;
    'Prediction ratio', @f2p2;
    'Average threshold', @f2p3;
    };
    
    'Others', {
    'Interleaved vs. combined alone', @f9p1;
    'Time-shifting analysis', @f9p2;
    'With or without tolerance', @f9p3;
    'Choice patterns', @f9p4;
    }
    
    };

%% ====================================== Function Definitions =============================================%

    function f1p1(debug)       % Threshold and Bias: Training day
        if debug;  dbstack;  keyboard;  end
        
        set(figure(76),'Pos',[10 258 880 705]); clf; hold on;
        
        if monkey == 2 % ~isempty(strfind(mat_address,'Messi'))
            plot(1:length(Polo),Polo(:,1),'--b','LineW',2);
            plot(1:length(Polo),Polo(:,3),'--r','LineW',2);
            plot(1:length(Polo),Polo(:,5),'--g','LineW',2);
            
            set(legend('Polo'));
            set(gca,'yscale','log','ytick',[1 2 3 4 5 10 100]);
            %     ,'yticklabel',[1 5 10 100]);
        end
        
        % Plot threshold
        for condition = 1:3
            h(condition) = plot(1:nDay,Thresh_day(:,condition),[colors{condition} 'o'],'markerfacecolor',colors{condition});
            
            if ~all(isnan(Thresh_day(:,condition)))
                xx = find(~isnan(Thresh_day(:,condition)),1,'first'):find(~isnan(Thresh_day(:,condition)),1,'last');
                yy = smooth(Thresh_day(:,condition),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
            
            %      axis([0 50 2 190]);
        end
        
        % Plot coherence
        h(4) = plot(1:nDay,Coh_day,'rs');
        
        legend(h,{'Vestibular','Visual','Combined','Coherence'});
        set(gca,'yscale','log','ytick',[1 5 10]);
        
        xlabel('Training Day');
        ylabel('Threshold'); xlim([0 nDay+2]);
        gg = gca;
        GrayGrid(gg);
        Annotation(1:nDay,tf_day,glass_day);
        SetFigure(20);

        %% -- Abs(bias) --
        set(figure(77),'Pos',[67 182 887 719]); clf; hold on;
        
        if monkey == 2 % ~isempty(strfind(mat_address,'Messi'))
            plot(1:length(Polo),abs(Polo(:,2)),'--b','LineW',2);
            set(legend('Polo'));
            set(gca,'yscale','log','ytick',[0.1 1 10],'yticklabel',[0.1 1 10]);
            ylim([0.01 40]);
        else
            set(gca,'yscale','log','ytick',[0.1 1 10],'yticklabel',[0.1 1 10]);
        end
        
        for condition = 1:3
            plot(1:nDay,abs(Bias_day(:,condition)),[colors{condition} 's'],'markerfacecolor',colors{condition});
            
            if ~all(isnan(Bias_day(:,condition)))
                xx = find(~isnan(Bias_day(:,condition)),1,'first'):find(~isnan(Bias_day(:,condition)),1,'last');
                yy = smooth(abs(Bias_day(:,condition)),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
        end
        
        xlim([0 nDay+2]);
        xlabel('Training Day');
        ylabel('|Bias|');
        
        GrayGrid(gca);
        
        SetFigure(20);
        
        %% == u time course ==
        
        set(figure(78),'position',[119 128 880 705]); clf
        plot([0 nDay],[0 0],'k--','linew',2);
        hold on
        
        for condition = 1:3
            plot(Bias_day(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10);
            
            if ~all(isnan(Bias_day(:,condition)))
                xx = find(~isnan(Bias_day(:,condition)),1,'first'):find(~isnan(Bias_day(:,condition)),1,'last');
                yy = smooth(Bias_day(:,condition),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
        end
        % plot(Bias_pred_session,'go','markersize',13,'linewidt',1.5);
        xlabel('Training Day');
        ylabel('Bias');
        xlim([0 nDay+2]);
        ylim(1.1*[-max(abs(ylim)) max(abs(ylim))]);
        
        Annotation(1:nDay,tf_day,glass_day);
        
        SetFigure();
    end
    function f1p2(debug)       % Prediction ratio : Training Day
        if debug;  dbstack;  keyboard;  end
        
        sigma_pred_day = 1./sqrt(1./Thresh_day(:,1).^2 + 1./Thresh_day(:,2).^2);
        sigma_pred_ratio_day = Thresh_day(:,3)./sigma_pred_day;
        
        set(figure(79),'position',[19 85 1388 535]);
        clf;
        plot([0 nDay],[1 1],'k--','linew',2); hold on;
        
        % Min threshold
        % plot(min(Thresh_session(:,1:2),[],2)./sigma_pred,'v','color',[0.5 0.5 0.5],'markerfacecol',[0.5 0.5 0.5]);
        % plot(max(Thresh_session(:,1:2),[],2)./sigma_pred,'k^','markerfacecol','k');
        
        if ~all(isnan(sigma_pred_day))
            plot(smooth(min(Thresh_day(:,1:2),[],2)./sigma_pred_day,20,'rloess'),['-.k'],'linewid',2);
            plot(smooth(max(Thresh_day(:,1:2),[],2)./sigma_pred_day,20,'rloess'),['-.k'],'linewid',2);
        end
        
        
        % plot(1:failureBegin-1,sigma_pred_ratio(1:failureBegin-1),'ks','markerfacecol','k','markersize',9);
        % plot(failureBegin:n, sigma_pred_ratio(failureBegin:end),'ks','markersize',9);
        
        if ~showPlatform
            scatter(1:length(sigma_pred_ratio_day),sigma_pred_ratio_day,150,linspace(0,1,length(sigma_pred_ratio_day)),'fill','s');
        else
            for pp = 1:length(platforms)
                ind = find(platform_day == platforms(pp));
                plot(ind,sigma_pred_ratio_day(ind),'s','color',platformColor{pp},'markerfacecol',platformColor{pp},'markersize',9);
            end
        end
        
        if ~all(isnan(sigma_pred_ratio_day))
            plot(smooth(sigma_pred_ratio_day,30,'rloess'),['-k'],'linewid',2);
        end
        
        ylim([0.3 3]); xlim([0 nDay+2]);
        set(gca,'YScale','log','Ytick',[0.5 1 2:3]);
        
        xlabel('Training day');
        ylabel('Prediction ratio');
        % ylabel('Prediction ratio = actual / predicted threshold');
        
        Annotation(1:nDay,tf_day,glass_day);
        
        SetFigure();
    end

    function f1p3(debug)       %  Averaged threshold
        if debug;  dbstack;  keyboard;  end
       
        figure(20); clf; hold on
        
        if monkey == 2 % Messi
            average_duration = 55:nDay;
        elseif monkey == 1 % Polo
            average_duration = 110:nDay;
        end
        
        data_to_average = Thresh_day(average_duration,:);
        data_to_average(any(isnan(data_to_average),2),:) = [];
        
        % Predicted threshold
        data_to_average(:,4) = sqrt(data_to_average(:,1).^2 .* data_to_average(:,2).^2 ./ (data_to_average(:,1).^2 + data_to_average(:,2).^2));
        
        mean_threshold = mean(data_to_average);
        sem_threshold = std(data_to_average)/sqrt(size(data_to_average,1));
        
        
        % Plotting
        
        color_bar = {'b','r','g','c'};
        
        xlim([0.5 4.5]);
        
        for i = 1:4
            bar(i,mean_threshold(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_threshold(i),sem_threshold(i),color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_vis] = ttest(data_to_average(:,1)-data_to_average(:,2));
        [~,p_comb_vest] = ttest(data_to_average(:,1)-data_to_average(:,3));
        [~,p_comb_vis] = ttest(data_to_average(:,2)-data_to_average(:,3));
        [~,p_comb_pred] = ttest(data_to_average(:,3)-data_to_average(:,4));
        
        set(gca,'xtick',1:4,'xticklabel',{'Vestibular','Visual','Combined','Optimal'});
        rotateXLabels(gca,45);
        ylabel('Threshold');
        
        text(0,max(ylim),sprintf('p vest-vis = %g, p comb-vest = %g\np comb-vis = %g, p comb-pred = %g',p_vest_vis,p_comb_vest,p_comb_vis,p_comb_pred),'fontsize',10.5);
        title(sprintf('n = %g days',length(data_to_average)));
        SetFigure(22);
    end

    function f2p1(debug)       %  Threshold and Bias: Session
        if debug;  dbstack;  keyboard;  end
        
        % == sigma time course ==
        
        set(figure(3),'position',[19 85 1388 535]); clf
        
        comb_alone =  isnan(Thresh_session(:,1)) & isnan(Thresh_session(:,2)) & ~isnan(Thresh_session(:,3)); % Combine alone vs interleaved
        
        for condition = 1:3
            plot(Thresh_session(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10); hold on;
            if ~all(isnan(Thresh_session(:,condition)))
                plot(smooth(Thresh_session(:,condition),30,'rloess'),['-' colors{condition}],'linewid',2);
            end
        end
        
        % Replot comb_alone
        plot(find(comb_alone),Thresh_session(comb_alone,3),'ko','markerfacecol','g','markersize',10,'linewi',2);
        
        % plot(sigma_pred_session,'go','markersize',10,'linewidt',1.5);
        
        % plot(coh);
        xlabel('Session');
        ylabel('Threshold');
        xlim([0 nSession+2]);
        set(gca,'YScale','log','Ytick',[1:10]);
        
        Annotation(day,tf_session,glass_session);
        
        SetFigure();
        
        % == u time course ==
        
        Bias_pred_session = (Thresh_session(:,2).^2 .* Bias_session(:,1) + Thresh_session(:,1).^2 .* Bias_session(:,2))./(Thresh_session(:,1).^2 + Thresh_session(:,2).^2);
        
        set(figure(4),'position',[19 85 1388 535]); clf
        plot([0 nSession],[0 0],'k--','linew',2);
        hold on
        
        for condition = 1:3
            plot(Bias_session(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10);
            if ~all(isnan(Bias_session(:,condition)))
                plot(smooth(Bias_session(:,condition),30,'rloess'),['-' colors{condition}],'linewid',2);
            end
        end
        
        % Replot comb_alone
        plot(find(comb_alone),Bias_session(comb_alone,3),'ko','markerfacecol','g','markersize',10,'linewi',2);
        
        % plot(Bias_pred_session,'go','markersize',13,'linewidt',1.5);
        xlabel('Session');
        ylabel('Bias');
        ylim(1.1*[-max(abs(ylim)) max(abs(ylim))]);
        xlim([0 nSession+2]);
        
        Annotation(day,tf_session,glass_session);
        
        
        SetFigure();
    end
    function f2p2(debug)       % Prediction ratio: Session
        if debug;  dbstack;  keyboard;  end
        
        sigma_pred_session = 1./sqrt(1./Thresh_session (:,1).^2 + 1./Thresh_session(:,2).^2);
        sigma_pred_ratio_session = Thresh_session(:,3)./sigma_pred_session;
        
        set(figure(5),'position',[19 85 1388 535]);
        clf;
        plot([0 nSession],[1 1],'k--','linew',2); hold on;
        
        % Min threshold
        % plot(min(Thresh_session(:,1:2),[],2)./sigma_pred,'v','color',[0.5 0.5 0.5],'markerfacecol',[0.5 0.5 0.5]);
        % plot(max(Thresh_session(:,1:2),[],2)./sigma_pred,'k^','markerfacecol','k');
        
        if ~all(isnan(sigma_pred_session))
            plot(smooth(min(Thresh_session(:,1:2),[],2)./sigma_pred_session,20,'rloess'),['-.k'],'linewid',2);
            plot(smooth(max(Thresh_session(:,1:2),[],2)./sigma_pred_session,20,'rloess'),['-.k'],'linewid',2);
        end
        
        % plot(1:failureBegin-1,sigma_pred_ratio(1:failureBegin-1),'ks','markerfacecol','k','markersize',9);
        % plot(failureBegin:n, sigma_pred_ratio(failureBegin:end),'ks','markersize',9);
        
        if ~showPlatform
            scatter(1:nSession,sigma_pred_ratio_session,150,linspace(0,1,nSession),'fill','s');
        else
            for pp = 1:length(platforms)
                ind = find(platform_session == platforms(pp));
                plot(ind,sigma_pred_ratio_session(ind),'s','color',platformColor{pp},'markerfacecol',platformColor{pp},'markersize',9);
            end
        end
        
        if ~all(isnan(sigma_pred_ratio_session))
            plot(smooth(sigma_pred_ratio_session,30,'rloess'),['-k'],'linewid',2);
        end
        
        ylim([0.3 3]); xlim([0 nSession+2]);
        set(gca,'YScale','log','Ytick',[0.5 1 2:3]);
        
        xlabel('Session');
        ylabel('Prediction ratio');
        % ylabel('Prediction ratio = actual / predicted threshold');
        
        Annotation(day,tf_session,glass_session);
        
        SetFigure();
    end

    function f2p3(debug)       %  Averaged threshold
        if debug;  dbstack;  keyboard;  end
        
        figure(20); clf; hold on
        
        if monkey == 2 % Messi
            average_duration = 180:nSession;
        elseif monkey == 1 % Polo
            average_duration = 160:nSession;
        end
        
        data_to_average = Thresh_session(average_duration,:);
        data_to_average(any(isnan(data_to_average),2),:) = [];
        
        % Predicted threshold
        data_to_average(:,4) = sqrt(data_to_average(:,1).^2 .* data_to_average(:,2).^2 ./ (data_to_average(:,1).^2 + data_to_average(:,2).^2));
        
        mean_threshold = mean(data_to_average);
        sem_threshold = std(data_to_average)/sqrt(size(data_to_average,1));
        
        
        % Plotting
        
        color_bar = {'b','r','g','c'};
        
        xlim([0.5 4.5]);
        
        for i = 1:4
            bar(i,mean_threshold(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_threshold(i),sem_threshold(i),color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_vis] = ttest(data_to_average(:,1)-data_to_average(:,2));
        [~,p_comb_vest] = ttest(data_to_average(:,1)-data_to_average(:,3));
        [~,p_comb_vis] = ttest(data_to_average(:,2)-data_to_average(:,3));
        [~,p_comb_pred] = ttest(data_to_average(:,3)-data_to_average(:,4));
        
        set(gca,'xtick',1:4,'xticklabel',{'Vestibular','Visual','Combined','Optimal'});
        rotateXLabels(gca,45);
        ylabel('Threshold');
        
        text(0,max(ylim),sprintf('p vest-vis = %g, p comb-vest = %g\np comb-vis = %g, p comb-pred = %g',p_vest_vis,p_comb_vest,p_comb_vis,p_comb_pred),'fontsize',10.5);
        title(sprintf('n = %g sessions',length(data_to_average)));
        SetFigure(22);
    end

    function f9p1(debug)       % Interaction analysis: interleaved vs comb_alone
        if debug;  dbstack;  keyboard;  end

        figure(10); clf; hold on
        
        means = mean(interaction_pairs,1);
        se = std(interaction_pairs,0,1)/sqrt(size(interaction_pairs,1));
        
        bar([1 2],means,0.5,'facecol','k','edgecol','none');
        h = errorbar([1 2],means,se,'k','linestyle','none','linewidth',5);
        errorbar_tick(h,5);
        
        plot(interaction_pairs','o-','color',[0.6 0.6 0.6],'markersize',7,'markerfacecol',[0.6 0.6 0.6],'linew',2);
        
        [~,p]=ttest(interaction_pairs(:,1),interaction_pairs(:,2));
        
        set(gca,'xtick',[1 2],'xticklabel',{'Comb alone','Interleaved'});
        text(1.5,1,sprintf('paired t-test\n \\itp\\rm = %3.3g',p));
        xlim([0.5 2.5]);
        ylabel('Actual / Predicted');
        
        SetFigure(20);
    end
    function f9p2(debug)       % Time-shifting-related
        if debug;  dbstack;  keyboard;  end

        % figure(88); clf; hold on;
        stim_type = {'Vest','Visual','Combined','Prediction ratio'};
        
        shift_ratio_all = cell(1,4);
        shift_ratio_all(:) = {nan(length(Thresh_shift_session),200)};
        
        z_score_time_shift_all = shift_ratio_all;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thresholdMask = 100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Prediction ratio shift
        for i = 1:nSession
            if ~any(isnan([Thresh_shift_session{i,:}]))  % Have all three conditions
                pred =  1./sqrt(1./Thresh_shift_session{i,1}.^2 + 1./Thresh_shift_session{i,2}.^2);
                Thresh_shift_session{i,4} = Thresh_shift_session{i,3}./ pred;
            else
                Thresh_shift_session{i,4} = nan;
            end
        end
        
        for condition = 1:4  % The last for prediction ratio
            for i = 1:nSession
                
                if sum(Thresh_shift_session{i,condition} < thresholdMask)>0
                    
                    Thresh_shift_session_masked = Thresh_shift_session{i,condition}(Thresh_shift_session{i,condition}<thresholdMask);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    shift_ratio = Thresh_shift_session_masked;
                    %             shift_ratio = Thresh_shift_session_masked/Thresh_shift_session_masked(1);
                    %             shift_ratio = Thresh_shift_session_masked/geomean(Thresh_shift_session_masked);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    shift_ratio_all{condition}(i,1:length(shift_ratio))=shift_ratio;
                    
                    %             plot(shift_ratio);
                    
                    % Z-score
                    mean_ = mean(Thresh_shift_session_masked);
                    std_ = std(Thresh_shift_session_masked);
                    z_score_time_shift = (Thresh_shift_session_masked - mean_)/std_ ;
                    z_score_time_shift_all{condition}(i,1:length(z_score_time_shift)) = z_score_time_shift;
                end
                
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nMask = 5;
            seg = 5; % Num of stages across all training session
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            maxTimeShift{condition} = find(sum(~isnan(shift_ratio_all{condition}))>nMask,1,'last');
            
            shift_ratio_mean = []; shift_ratio_se = []; shift_ratio_mean_seg = [];
            
            for j = 1:maxTimeShift{condition}
                xx = shift_ratio_all{condition}(~isnan(shift_ratio_all{condition}(:,j)),j);
                shift_ratio_mean(j) = mean(log(xx));
                shift_ratio_se(j) = std(log(xx))/sqrt(length(xx));
                
                for ss = 1:seg
                    shift_ratio_all_this_seg = shift_ratio_all{condition}(1+fix((ss-1)/seg*end):fix(ss/seg*end),j);
                    shift_ratio_mean_seg(j,ss) = mean(log(shift_ratio_all_this_seg(~isnan(shift_ratio_all_this_seg))));
                end
            end
            
            figure(89+condition); clf;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     set(gca,'ColorOrder',1-0.9./repmat([seg:-1:1]',1,3));
            co = colormap(spring);
            set(gca,'ColorOrder',co(round(linspace(1,size(co,1),seg)),:));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            hold on;
            plot(exp(shift_ratio_mean_seg),'LineWidth',3);
            errorbar(1:maxTimeShift{condition},exp(shift_ratio_mean),exp(shift_ratio_mean)-exp(shift_ratio_mean-shift_ratio_se),exp(shift_ratio_mean+shift_ratio_se)-exp(shift_ratio_mean),'k');
            %     axis([0 50 0.7 1.1]);
            plot(xlim,[1 1],'k--');
            xlabel('Repetition'); ylabel('Normalized threshold');
            SetFigure(20);
            title(stim_type{condition});
        end
        
    end
    function f9p3(debug)       % With or without tolerance
        if debug;  dbstack;  keyboard;  end

        Thresh_session_no_tol = reshape(cell2mat([group_result.Thresh_psy]),3,[])';
        Thresh_shift_session_no_tol = reshape([group_result.psy_thresh_shift],3,[])';
        
        Thresh_session_tol = reshape(cell2mat([group_result.Thresh_psy_tol]),3,[])';
        Thresh_shift_session_tol = reshape([group_result.psy_thresh_shift_tol],3,[])';
                
        sigma_pred_session_no_tol = 1./sqrt(1./Thresh_session_no_tol (:,1).^2 + 1./Thresh_session_no_tol(:,2).^2);
        sigma_pred_ratio_session_no_tol = Thresh_session_no_tol(:,3)./sigma_pred_session_no_tol;
        
        sigma_pred_session_tol = 1./sqrt(1./Thresh_session_tol (:,1).^2 + 1./Thresh_session_tol(:,2).^2);
        sigma_pred_ratio_session_tol = Thresh_session_tol(:,3)./sigma_pred_session_tol;
        
        figure(930); clf;  plot(1:nSession, sigma_pred_ratio_session_tol,'r-o'); hold on;
        plot(1:nSession,sigma_pred_ratio_session_no_tol','-o');  axis tight;
        legend('tol','No tol');        SetFigure(20);
        
        %% Comparison
        sigma_pred_ratio_pairs = [sigma_pred_ratio_session_no_tol sigma_pred_ratio_session_tol];
        sigma_pred_ratio_pairs(abs(diff(sigma_pred_ratio_pairs,[],2)) < eps,:) = nan;
        
        figure(931); clf; plot(sigma_pred_ratio_pairs(:,1),sigma_pred_ratio_pairs(:,2),'ro');  axis square; hold on;
        mins = min([xlim ylim]); maxs = max([xlim ylim]);  axis([mins maxs mins maxs]);
        plot([mins maxs], [mins maxs],'k--'); xlabel('Pred ratio no tol'); ylabel('Pred ratio tol'); SetFigure(20);
        
        figure(932); clf; hold on
        
        means = nanmean(sigma_pred_ratio_pairs,1);
        se = nanstd(sigma_pred_ratio_pairs,0,1)/sqrt(sum(~isnan(sigma_pred_ratio_pairs(:,1)),1));
        
        bar([1 2],means,0.5,'facecol','k','edgecol','none');  hold on; 
        h = errorbar([1 2],means,se,'k','linestyle','none','linewidth',5); 
        errorbar_tick(h,5);
        
        plot(sigma_pred_ratio_pairs','o-','color',[0.6 0.6 0.6],'markersize',7,'markerfacecol',[0.6 0.6 0.6],'linew',2);
        
        [~,p]=ttest(sigma_pred_ratio_pairs(:,1),sigma_pred_ratio_pairs(:,2));
        
        set(gca,'xtick',[1 2],'xticklabel',{'No Tol','Tol'});
        text(max(xlim)/2,max(ylim)*0.9,sprintf('paired t-test\n \\itp\\rm = %3.3g',p));
        xlim([0.5 2.5]);
        ylabel('Prediction ratio');
        
        SetFigure(20);

    end
    function f9p4(debug)       % Choice patterns
        if debug;  dbstack;  keyboard;  end

        for i=1:nSession
            raw = group_result(i).raw;
            
            % === All angles ===
            win = raw(1:end-1,5) == 0;
            lose = raw(1:end-1,5) == 5;
            stay = raw(2:end,4) == raw(1:end-1,4);
            shift = raw(2:end,4) ~= raw(1:end-1,4);
            
            win_shift(i,:) = sum(win & shift)/sum(win);
            p_win_shift(i,:) = p_value_binomial(sum(win),sum(win & shift));
            lose_shift(i,:) = sum(lose & shift)/sum(lose);
            p_lose_shift(i,:) = p_value_binomial(sum(lose),sum(lose & shift));
            
            all_shift(i,:) = sum(shift)/(length(raw)-1);
            p_all_shift(i,:) = p_value_binomial(length(raw)-1,sum(shift));
            
            % Heading shift that it actually was
            heading_shift(i,:) = sum(sign(raw(1:end-1,3)) ~= sign(raw(2:end,3)))/(size(raw,1)-1);
            
            % ======  Shift prop V.S. CR =====
            
            % CR of different abs(heading)
            unique_abs_heading{i} = unique(abs(raw(:,3)));
            
            % Shift prop per abs(heading)
            for ah = 1:length(unique_abs_heading{i})
                this_heading = find(abs(raw(:,3)) == unique_abs_heading{i}(ah));
                this_heading(this_heading==1) = []; % Ignore the first trial
                
                win_abs_heading_last = raw(this_heading-1,5) == 0;
                lose_abs_heading_last = raw(this_heading-1,5) == 5;
                
                win_abs_heading_this = raw(this_heading,5) == 0;
                lose_abs_heading_this = raw(this_heading,5) == 5;
                
                shift_abs_heading_last = raw(this_heading-1,4) ~= raw(this_heading,4);
                
                CR_abs_heading{i}(ah) = sum(win_abs_heading_this) / length(this_heading);
                
                win_shift_abs_heading_last{i}(ah) = sum(win_abs_heading_last & shift_abs_heading_last)/sum(win_abs_heading_last);
                p_win_shift_abs_heading_last{i}(ah) = p_value_binomial(sum(win_abs_heading_last),sum(win_abs_heading_last & shift_abs_heading_last),0.5);
                
                if sum(lose_abs_heading_last)<1
                    lose_shift_abs_heading_last{i}(ah) = nan;
                    p_lose_shift_abs_heading_last{i}(ah) = nan;
                else
                    lose_shift_abs_heading_last{i}(ah) = sum(lose_abs_heading_last & shift_abs_heading_last)/sum(lose_abs_heading_last);
                    p_lose_shift_abs_heading_last{i}(ah) = p_value_binomial(sum(lose_abs_heading_last),sum(lose_abs_heading_last & shift_abs_heading_last),0.5 );
                end
                
                win_shift_abs_heading_this{i}(ah) = sum(win_abs_heading_this & shift_abs_heading_last)/sum(win_abs_heading_this);
                p_win_shift_abs_heading_this{i}(ah) = p_value_binomial(sum(win_abs_heading_this),sum(win_abs_heading_this & shift_abs_heading_last),0.5);
                
                if sum(lose_abs_heading_this)<1
                    lose_shift_abs_heading_this{i}(ah) = nan;
                    p_lose_shift_abs_heading_this{i}(ah) = nan;
                else
                    lose_shift_abs_heading_this{i}(ah) = sum(lose_abs_heading_this & shift_abs_heading_last)/sum(lose_abs_heading_this);
                    p_lose_shift_abs_heading_this{i}(ah) = p_value_binomial(sum(lose_abs_heading_this),sum(lose_abs_heading_this & shift_abs_heading_last),0.5);
                end
                
            end
            
            %     % ====== Using HD constant per se (permutation test) ======
            %     if ~staircase(i) % Only apply to constant sessions
            %         output = p_value_constantHD(length(unique(group_result(i).fit_data_psycho_cum{1}(:,1))),...
            %             group_result(i).repetitionN,[win_shift(i,:) lose_shift(i,:) all_shift(i,:)]);
            %         constant_per_se(i,:) = output.mean;
            %         p_win_shift_constant_per_se(i,:) = output.p(1);
            %         p_lose_shift_constant_per_se(i,:) = output.p(2);
            %         p_all_shift_constant_per_se(i,:) = output.p(3);
            %     else
            %         constant_per_se(i,:) = nan;
            %         p_win_shift_constant_per_se(i,:) = nan;
            %         p_lose_shift_constant_per_se(i,:) = nan;
            %         p_all_shift_constant_per_se(i,:) = nan;
            %         p_all_constant_per_se(i,:) = nan;
            %     end
            
            
        end
        
        %% Time course
        
        both_sign = p_win_shift < 0.05 & p_lose_shift < 0.05;
        none_sign = p_win_shift >= 0.05 & p_lose_shift >= 0.05;
        one_sign = xor(p_win_shift<0.05,p_lose_shift<0.05);
        
        set(figure(78),'Pos',[18 84 838 458]);  clf; subplot(2,1,1);
        % plot(win_shift,'k');
        hold on;
        plot(find(p_win_shift<0.05 & staircase),win_shift(p_win_shift<0.05 & staircase),'ro','markerfacecolor','r');
        plot(find(p_win_shift>=0.05& staircase),win_shift(p_win_shift>=0.05& staircase),'ro','markerfacecolor','none');
        plot(find(p_win_shift<0.05& ~staircase),win_shift(p_win_shift<0.05& ~staircase),'ko','markerfacecolor','k');
        plot(find(p_win_shift>=0.05& ~staircase),win_shift(p_win_shift>=0.05& ~staircase),'ko','markerfacecolor','none');
        
        plot(xlim,[0.5 0.5],'k--'); ylabel('Win-shift'); axis tight; ylim(0.5+1.1*[-max(abs(ylim-0.5)) max(abs(ylim-0.5))]);
        
        subplot(2,1,2);
        % plot(lose_shift,'k');
        hold on;
        plot(find(p_lose_shift<0.05 & staircase),lose_shift(p_lose_shift<0.05 & staircase),'ro','markerfacecolor','r');
        plot(find(p_lose_shift>=0.05& staircase),lose_shift(p_lose_shift>=0.05& staircase),'ro','markerfacecolor','none');
        plot(find(p_lose_shift<0.05& ~staircase),lose_shift(p_lose_shift<0.05& ~staircase),'ko','markerfacecolor','k');
        plot(find(p_lose_shift>=0.05& ~staircase),lose_shift(p_lose_shift>=0.05& ~staircase),'ko','markerfacecolor','none');
        plot(xlim,[0.5 0.5],'k--'); ylabel('Lose-shift');   axis tight;  ylim(0.5+1.1*[-max(abs(ylim-0.5)) max(abs(ylim-0.5))]);
        SetFigure(20);
        set(findall(gcf,'type','line'),'markersize',10)
        set(findall(gcf,'type','line'),'linewi',1.5)
        
        %% Histograms
        set(figure(2699),'Posi',[708 189 867 705]); clf;
        xCenters = [0.2:0.02:0.8];
        
        hh(1)= HistComparison({win_shift(p_win_shift<0.05 & staircase),win_shift(p_win_shift>=0.05 & staircase)},...
            'TTest',0.5,'EdgeColors',{'r','r'},'FaceColors',{'r','none'},'XCenters',xCenters,'Axes',subplot(221)); title('Win shift staircase');
        hh(2)= HistComparison({lose_shift(p_lose_shift<0.05 & staircase),lose_shift(p_lose_shift>=0.05 & staircase)},...
            'TTest',0.5,'EdgeColors',{'r','r'},'FaceColors',{'r','none'},'XCenters',xCenters,'Axes',subplot(223)); title('Lose shift staircase');
        hh(3)= HistComparison({win_shift(p_win_shift<0.05 & ~staircase),win_shift(p_win_shift>=0.05 & ~staircase)},...
            'TTest',0.5,'XCenters',xCenters,'Axes',subplot(222)); title('Win shift constant');
        hh(4)= HistComparison({lose_shift(p_lose_shift<0.05 & ~staircase),lose_shift(p_lose_shift>=0.05 & ~staircase)},...
            'TTest',0.5,'XCenters',xCenters,'Axes',subplot(224)); title('Lose shift constant');
        
        set([hh.ax],'Xlim',[0.2 0.8])
        
        
        % % HD Constant shift prop per se
        % load('HDConstantShiftProp.mat');
        %
        % h = HistComparison({heading_shift(staircase),heading_shift(~staircase)},'Style','grouped',...
        %     'TTest',0.5,'EdgeColors',{'r','k'},'FaceColors',{'r','k'},'XNumbers',20); xlim([0.2 0.8]);
        % hold on;
        %
        % % xx = 0:0.005:1;
        % % plot(xx,fitresults{4}(xx)*length(heading_shift(~staircase))*(h.x_centers(2)-h.x_centers(1)),'linewidth',2);
        % %
        % title('Heading shift');
        % SetFigure(20);
        
        %% ---- CR v.s. Shift -----
        
        sig = [p_lose_shift_abs_heading_last{~staircase}]'<0.05; nsig = [p_lose_shift_abs_heading_last{~staircase}]'>=0.05;
        CRs = [CR_abs_heading{~staircase}]'; props = [lose_shift_abs_heading_last{~staircase}]';
        h=LinearCorrelation({CRs(sig),CRs(nsig)},{props(sig),props(nsig)},'Method','Pearson','FittingMethod',2,...
            'FigN',179,'XLabel','Correct rate','YLabel','P ( Shift | Last-lose)','YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','CombinedIndex',[3]);
        plot(xlim,[0.5 0.5],'k--'); ylim([0 1]); set(h.ax_yhist,'xlim',ylim); %delete([h.group(1:2).line]);
        
        sig = [p_win_shift_abs_heading_last{~staircase}]'<0.05; nsig = [p_win_shift_abs_heading_last{~staircase}]'>=0.05;
        CRs = [CR_abs_heading{~staircase}]'; props = [win_shift_abs_heading_last{~staircase}]';
        h=LinearCorrelation({CRs(sig),CRs(nsig)},{props(sig),props(nsig)},'Method','Pearson','FittingMethod',2,...
            'FigN',180,'XLabel','Correct rate','YLabel','P ( Shift | Last-win)','YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','CombinedIndex',[3]);
        plot(xlim,[0.5 0.5],'k--'); ylim([0 1]); set(h.ax_yhist,'xlim',ylim); %delete([h.group(1:2).line]);
        
        sig = [p_lose_shift_abs_heading_this{~staircase}]'<0.05; nsig = [p_lose_shift_abs_heading_this{~staircase}]'>=0.05;
        CRs = [CR_abs_heading{~staircase}]'; props = [lose_shift_abs_heading_this{~staircase}]';
        h=LinearCorrelation({CRs(sig),CRs(nsig)},{props(sig),props(nsig)},'Method','Pearson','FittingMethod',2,...
            'FigN',181,'XLabel','Correct rate','YLabel','P ( Shift | This-lose)','YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','CombinedIndex',[3]);
        plot(xlim,[0.5 0.5],'k--'); ylim([0 1]); set(h.ax_yhist,'xlim',ylim);  %delete([h.group(1:2).line]);
        
        sig = [p_win_shift_abs_heading_this{~staircase}]'<0.05; nsig = [p_win_shift_abs_heading_this{~staircase}]'>=0.05;
        CRs = [CR_abs_heading{~staircase}]'; props = [win_shift_abs_heading_this{~staircase}]';
        h=LinearCorrelation({CRs(sig),CRs(nsig)},{props(sig),props(nsig)},'Method','Pearson','FittingMethod',2,...
            'FigN',182,'XLabel','Correct rate','YLabel','P ( Shift | This-win)','YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','CombinedIndex',[3]);
        plot(xlim,[0.5 0.5],'k--'); ylim([0 1]); set(h.ax_yhist,'xlim',ylim); %delete([h.group(1:2).line]);
        
        % % Per se (win shift)
        % set(figure(79),'Pos',[18 84 838 458]);  clf; subplot(2,1,1);
        % plot(win_shift,'k'); hold on;
        % plot(constant_per_se,'sb');
        % plot(find(p_win_shift_constant_per_se<0.05 & ~staircase),win_shift(p_win_shift_constant_per_se<0.05 & ~staircase),'ko','markerfacecolor','k');
        % plot(find(p_win_shift_constant_per_se>=0.05 & ~staircase),win_shift(p_win_shift_constant_per_se>=0.05 & ~staircase),'ko','markerfacecolor','none');
        % plot(xlim,[0.5 0.5],'k--'); ylabel('Win Shift');
        % HistComparison({win_shift(p_win_shift_constant_per_se<0.05& ~staircase),win_shift(p_win_shift_constant_per_se>=0.05& ~staircase)},...
        %     'TTest',mean(constant_per_se(~staircase)),'XCenters',[0.2:0.02:0.8],'Axes',subplot(223)); title('Win shift constant (Per se)');
        % SetFigure(15);
        
        
        % Plotting
        h=LinearCorrelation({win_shift(none_sign& staircase),win_shift(one_sign& staircase),win_shift(both_sign& staircase),...
            win_shift(none_sign& ~staircase),win_shift(one_sign& ~staircase),win_shift(both_sign& ~staircase)},...
            {lose_shift(none_sign& staircase),lose_shift(one_sign& staircase),lose_shift(both_sign& staircase),...
            lose_shift(none_sign& ~staircase),lose_shift(one_sign& ~staircase),lose_shift(both_sign& ~staircase)},...
            'LineStyles',{'r:','r:','r:','k:','k:','k:','r:','k:','k-'},'Markers',{'o'},'FaceColors',{'none',[1 0.8 1],'r','none',[.8 .8 .8],'k'},...
            'CombinedIndex',[7 56 63],'xlabel','Win Shift','ylabel','Lose Shift','SameScale',1,...
            'XHist',20,'YHist',30,'XHistStyle','stacked','YHistStyle','stacked','FigN',2501);
        
        delete([h.group(1:6).line]);
        axes(h.ax_raw);
        plot(xlim,[0.5 0.5],'k--');
        plot([0.5 0.5],ylim,'k--');
        
    end


%% Improvement v.s. correct rate
% correct_rate_t = CR_session(Thresh_session<100);
% Thresh_psy_t = Thresh_session(Thresh_session<100);
%
% h=LinearCorrelation(correct_rate_t(1:end-1),Thresh_psy_t(1:end-1)./Thresh_psy_t(2:end),...
%     'XHist',20,'YHist',20,'MarkerSize',10,'XLabel','Correct rate','YLabel','Next improvement','logy',1,'Method','Pearson','FittingMethod',2,'FigN',2499);
% axes(h.ax_raw); plot(xlim,[0 0],'k--');
%
% h=LinearCorrelation(Thresh_psy_t(1:end-1),Thresh_psy_t(1:end-1)./Thresh_psy_t(2:end),...
%     'XHist',20,'YHist',20,'MarkerSize',10,'XLabel','Threshold','YLabel','Next improvement','logy',1,'FigN',2500,'Method','Pearson','FittingMethod',2);
% axes(h.ax_raw); plot(xlim,[0 0],'k--');


    function p = p_value_binomial(totalN, n, prop)
        if nargin < 3
            prop = 0.5;
        end
        
        x = 0:1:totalN;
        y = binopdf(x,totalN,prop);
        
        if n < 1  % Prop
            n = totalN * n;
        end
        
        % Putting this ahead makes the p_value symmetric (espeically when totalN is small).  HH20150125
        if n > totalN/2
            n = totalN-n;
        end
        
        p = sum(y(1:min(find(x>=n))));
        
    end

    function output = p_value_constantHD(headingN, repN, prop)
        simN = 1500;
        xcenterN = 200;
        
        unique_heading = [-fix(headingN/2):-1 , 1:fix(headingN/2)];
        
        headings = zeros(length(unique_heading),repN*simN);
        
        for i = 1:repN*simN
            headings(:,i) = unique_heading(randperm(length(unique_heading)));
        end
        
        headings = reshape(headings,length(unique_heading)*repN,[]);
        shift_prop = sum(sign(headings(1:end-1,:)) ~= sign(headings(2:end,:)))/(size(headings,1)-1);
        
        [n,x] = hist(shift_prop,xcenterN);
        
        output.mean = mean(shift_prop);
        
        for pp = 1:length(prop)
            output.p(pp) = sum(n(1:min(find(x>prop(pp)))))/ simN;
            if prop >= output.mean
                output.p = 1-output.p;
            end
        end
        
    end

    function Annotation(date,tf,glass)
        
        % Same day
        diff0 = find(diff(date)== 0);
        
        if ~isempty(diff0)
            diff0_begs = [diff0(1); diff0(1+find(diff(diff0)>1))];
            diff0_ends = [diff0(diff(diff0)>1) ;diff0(end)]+1;
        else
            diff0_begs = [];
        end
        
        ylims = ylim;
        xlims = xlim;
        
        for j = 1:length(diff0_begs);
            plot([diff0_begs(j)-0.2 diff0_ends(j)+0.2],[ylims(1) ylims(1)],'k-','linewid',10);
        end
        
        % Period that Polo refused to work
        % plot([70 70],ylims,'--k','linewi',2);
        % plot([74 74],ylims,'--k','linewi',2);
        
        
        % -- Target first v.s. 3D glass --
        a1 = gca;
        pos = get(gca,'Position');
        a2 = axes('position',[pos(1) pos(2)+pos(4) pos(3) 0.03]);
        hold on;
        
        if ~all(isnan(tf))
            plot(find(~isnan(tf) & tf ~= 0),0,'cs','markerfacecolor','c');
        end
        
        if ~all(isnan(glass)) && ~isempty(find(glass>0))
            plot(find(glass>0),1,'ms','markerfacecolor','m');
        end
        
        axis off;
        xlim(xlims); ylim([0 1])
        
        hlink = linkprop([a1 a2],'xlim');
        set(gcf,'UserData',hlink); % Store the hlink variable in an object's UserData property to keep the links! HH20150609
        
    end

end