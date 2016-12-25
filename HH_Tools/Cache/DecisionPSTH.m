% function DecisionPSTH
%%
[num,txt,raw] = xlsread('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',6);

%%
HEADS_N = 2;

AREA = 1;
REP = 3;
PSTH = 5:10;
CP = [22 35 48];
NT = CP-1;

psth = XlsCell2Mat(raw(HEADS_N+1:end,PSTH));
psth_t = str2num(raw{3,4});
cp_shift = XlsCell2Mat(raw(HEADS_N+1:end,CP));
nt_shift = XlsCell2Mat(raw(HEADS_N+1:end,NT));
cp_t = str2num(raw{end,20});
nt_t = cp_t;


areas = {'VIP','MST','LIP'};
areas_color = {'r','b','g'};
stim_type = {'Vestibular','Vis','Combined'};
stim_type_color = {'b','r','g'};

transparent = 0;

%%  Stim-type,  Different area

for k = 1:3
    
    % ------------------------- PSTH
    figure(400+k); clf; set(gcf,'position',[90 300 880 450]);  h_lines = [];
    for a = 1:3
        areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});
        pref = psth{(k-1)*2+1}(areaInd,:);
        null = psth{(k-1)*2+2}(areaInd,:);
        % raw
        %         plot(psth_t',psth{(k-1)*2+1}(areaInd,:)',[areas_color{a} '-']);  hold on; % Pref
        %         plot(psth_t',psth{(k-1)*2+2}(areaInd,:)',[areas_color{a} '--']); % Null
        
        % Normalized
        norm_range = find(psth_t>=-1800,1):find(psth_t>=0,1);
        
        offset = repmat(min(pref(:,norm_range),[],2),1,size(pref,2));
        gain = repmat(range(pref(:,norm_range),2),1,size(pref,2));
        norm_pref = (pref-offset)./gain; norm_pref(isnan(norm_pref(:,1)),:)=[];
        norm_null = (null-offset)./gain; norm_null(isnan(norm_null(:,1)),:)=[];
        
        %         plot(psth_t',norm_pref',[areas_color{a} '-']);  hold on; % Pref
        %         plot(psth_t',norm_null',[areas_color{a} '--']); % Null
        if ~isempty(norm_pref)
            h = shadedErrorBar(psth_t,norm_pref,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[areas_color{a} '-'],'linewidth',2},transparent);  hold on; % Pref
            h_lines(a) = h.mainLine;
            legends_lines{a} = [areas{a} ', N = ' num2str(size(norm_pref,1))];
        end
        if ~isempty(norm_null)
            shadedErrorBar(psth_t,norm_null,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[areas_color{a} '--'],'linewidth',2},transparent); % Null
        end
        axis tight
    end
    title([stim_type{k} '  PSTH']);
    ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
    legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
    xlabel('Time to Sac ON (ms)');
    ylabel('Normalized PSTH');
    SetFigure(15)
    
    % --------------------- CP
    figure(410+k); clf; set(gcf,'position',[90 300 880 450]);     h_lines = [];
    for a = 1:3
        areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});
        cp = cp_shift{k}(areaInd,:);
        % raw
        %         plot(psth_t',psth{(k-1)*2+1}(areaInd,:)',[areas_color{a} '-']);  hold on; % Pref
        %         plot(psth_t',psth{(k-1)*2+2}(areaInd,:)',[areas_color{a} '--']); % Null
        
        %         plot(cp',cp',[areas_color{a} '-']);  hold on; % Pref
        %         plot(cp_t',cp',[areas_color{a} '--']); % Null
        cp(isnan(cp(:,1)),:)=[];
        
        if ~isempty(cp)
            h = shadedErrorBar(cp_t,cp,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[areas_color{a} '-'],'linewidth',2},transparent);  hold on; % Pref
            h_lines(a) = h.mainLine;
            legends_lines{a} = [areas{a} ', N = ' num2str(size(cp,1))];
        end
        axis tight
    end
    title([stim_type{k} '   Grand CP']);
    xlabel('Center of 1000 ms window (to Stim ON)');
    ylabel('Grand CP');
    ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
    legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
    SetFigure(15)
    
    % --------------------- Neurothreshold
    figure(420+k); clf; set(gcf,'position',[90 300 880 450]);
    h_lines = [];
    for a = 1:3
        areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});
        nt = nt_shift{k}(areaInd,:);
        % raw
        %         plot(psth_t',psth{(k-1)*2+1}(areaInd,:)',[areas_color{a} '-']);  hold on; % Pref
        %         plot(psth_t',psth{(k-1)*2+2}(areaInd,:)',[areas_color{a} '--']); % Null
        
        %         plot(nt',nt',[areas_color{a} '-']);  hold on; % Pref
        %         plot(nt_t',nt',[areas_color{a} '--']); % Null
        nt(isnan(nt(:,1)),:)=[];
        
        if ~isempty(nt)
            nt_mean = exp(mean(log(nt),1));
            nt_upper = exp(mean(log(nt),1) + std(log(nt),0,1)/sqrt(size(nt,1)));
            nt_lower = exp(mean(log(nt),1) - std(log(nt),0,1)/sqrt(size(nt,1)));
            h = shadedErrorBar(cp_t,nt_mean, [nt_upper-nt_mean; nt_mean-nt_lower],{[areas_color{a} '-'],'linewidth',2},transparent);  hold on; % Pref
            %             plot(nt_t,nt_mean);hold on;
            h_lines(a) = h.mainLine;
            legends_lines{a} = [areas{a} ', N = ' num2str(size(nt,1))];
        end
        
    end
    axis tight
    ylim([0 80])
    xlabel('Center of 1000 ms window (to Stim ON)');
    ylabel('Neurothreshold');
    ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
    legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
    title([stim_type{k} '    Neurothreshold']);
    SetFigure(15)
    %     set(gca,'yscale','log','ytickmode','auto')
    
end


%% ------------------------- LIP: 3 stim-types

% ---------------------------------PSTH

figure(420); clf; set(gcf,'position',[90 300 880 450]);  h_lines = [];
a = 3;  % VIP only
areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});

all_psth = cell2mat(psth);
LIP_psth = all_psth(areaInd,:);

% Normalize psth according to visual response
norm_range = (find(psth_t>=-1800,1):find(psth_t>=0,1)) + (2-1)*2 * length(psth_t); % Visual
offset = repmat(min(LIP_psth(:,norm_range),[],2),1,size(LIP_psth,2));
gain = repmat(range(LIP_psth(:,norm_range),2),1,size(LIP_psth,2));
LIP_psth_norm = (LIP_psth-offset)./gain;

for k = 1:3
    norm_pref = LIP_psth_norm(:,(k-1)*2*length(psth_t)+1 : (k-1)*2*length(psth_t)+ length(psth_t));
    norm_null = LIP_psth_norm(:,((k-1)*2+1)*length(psth_t)+1 : ((k-1)*2+1)*length(psth_t)+ length(psth_t));
    % raw
%             plot(psth_t',norm_pref',[stim_type_color{k} '-']);  hold on; % Pref
%             plot(psth_t',norm_pref',[stim_type_color{k} '--']); % Null
%     
    norm_pref(isnan(norm_pref(:,1)),:)=[];
    norm_null(isnan(norm_null(:,1)),:)=[];
    
    %         plot(psth_t',norm_pref',[areas_color{a} '-']);  hold on; % Pref
    %         plot(psth_t',norm_null',[areas_color{a} '--']); % Null
    if ~isempty(norm_pref)
        h = shadedErrorBar(psth_t,norm_pref,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[stim_type_color{k} '-'],'linewidth',2},transparent);  hold on; % Pref
        h_lines(k) = h.mainLine;
        legends_lines{k} = [stim_type{k} ', N = ' num2str(size(norm_pref,1))];
    end
    if ~isempty(norm_null)
        shadedErrorBar(psth_t,norm_null,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[stim_type_color{k} '--'],'linewidth',2},transparent); % Null
    end
    axis tight
end
title('LIP, PSTH, different stim types');
ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
xlabel('Time to Sac ON (ms)');
ylabel('Normalized PSTH');
SetFigure(15)

% --------------------------------- CP

figure(421); clf; set(gcf,'position',[90 300 880 450]);  h_lines = [];
a = 3;  % VIP only
areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});

all_cp = cell2mat(cp_shift);
LIP_cp = all_cp(areaInd,:);

for k = 1:3
    norm_cp = LIP_cp(:,(k-1)*length(cp_t)+1 : (k-1)*length(cp_t)+ length(cp_t));
    % raw
%             plot(cp_t',norm_pref',[stim_type_color{k} '-']);  hold on; % Pref
%             plot(cp_t',norm_pref',[stim_type_color{k} '--']); % Null
%     
    
    %         plot(cp_t',norm_pref',[areas_color{a} '-']);  hold on; % Pref
    %         plot(cp_t',norm_null',[areas_color{a} '--']); % Null
    norm_cp(isnan(norm_cp(:,1)),:)=[];

    
    if ~isempty(norm_cp)
        h = shadedErrorBar(cp_t,norm_cp,{@(x) mean(x,1), @(x) std(x,0,1)/sqrt(size(x,1))},{[stim_type_color{k} '-'],'linewidth',2},transparent);  hold on; % Pref
        h_lines(k) = h.mainLine;
        legends_lines{k} = [stim_type{k} ', N = ' num2str(size(norm_cp,1))];
    end
    axis tight
end
title('LIP, Grand CP, different stim types');
ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
xlabel('Center of 1000 ms window (to Stim ON)');
ylabel('Grand CP');
SetFigure(15)

% --------------------------------- NT


figure(422); clf; set(gcf,'position',[90 300 880 450]);  h_lines = [];
a = 3;  % VIP only
areaInd = strcmp(txt(HEADS_N+1:end,1),areas{a});

all_nt = cell2mat(nt_shift);
LIP_nt = all_nt(areaInd,:);

for k = 1:3
    norm_nt = LIP_nt(:,(k-1)*length(nt_t)+1 : (k-1)*length(nt_t)+ length(nt_t));
    % raw
%             plot(nt_t',norm_pref',[stim_type_color{k} '-']);  hold on; % Pref
%             plot(nt_t',norm_pref',[stim_type_color{k} '--']); % Null
%     
    
    %         plot(nt_t',norm_pref',[areas_color{a} '-']);  hold on; % Pref
    %         plot(nt_t',norm_null',[areas_color{a} '--']); % Null
    norm_nt(isnan(norm_nt(:,1)),:)=[];
    
    if ~isempty(norm_nt)
        nt_mean = exp(mean(log(norm_nt),1));
        nt_upper = exp(mean(log(norm_nt),1) + std(log(norm_nt),0,1)/sqrt(size(norm_nt,1)));
        nt_lower = exp(mean(log(norm_nt),1) - std(log(norm_nt),0,1)/sqrt(size(norm_nt,1)));
        
        h = shadedErrorBar(nt_t,nt_mean, [nt_upper-nt_mean; nt_mean-nt_lower],{[stim_type_color{k} '-'],'linewidth',2},transparent);  hold on; % Pref
        
        h_lines(k) = h.mainLine;
        legends_lines{k} = [stim_type{k} ', N = ' num2str(size(norm_nt,1))];
    end
    axis tight
end
title('LIP, Neurothreshold, different stim types');
ylims = ylim; plot([0 0],ylims,'k-','linewidth',1.5)
legend(h_lines(h_lines~=0),legends_lines{h_lines~=0},'location','best');
xlabel('Center of 1000 ms window (to Stim ON)');
ylabel('Neurothreshold');
SetFigure(15)

% 
% 
% function out = normalize(raw,sel)
% offset = repmat(min(raw(:,sel),[],2),1,size(raw,2));
% gain = repmat(range(raw(:,sel),2),1,size(raw,2));
% out = (raw-offset)./gain;

