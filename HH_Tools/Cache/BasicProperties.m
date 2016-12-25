% function readXLS
%% Definition. Temporality hard coded
clear
headerWid = 3;
global xlsHeaders;
xlsHeaders = {  % Name, Loc in xls, length
    'GRIDLOC_', 3, 2;
    'SESSION_',3,1,
    'LOCX_', 5, 1;
    'LOCY_', 6, 1;
    'DEPTH_', 9, 1;
    'TYPE_', 10, 1;
    'FILE_',12,-1;
    'UNITS_', 21, 1;
    'Protocol_',23, 1;
    'SNR_', 27, 1;
    'SPKWAVE_', 28, -1;
    'RF_', 29, 4;
    'PREF_', 40, 1;
    'DDI_', 41, 1;
    'HTI_', 42, 1;
    'TUNING_P_', 43,1;
    'D90_',44,1;
    'D270_', 45, 1;
    'TUNING_MEAN_VIS_', 48, -1;
    'TUNING_MEAN_VEST_',38, -1;
    'TUNING_ERR_', 49, -1;
    
    % Del-Sac
    'LS_P_',67,1;
    'POSTS_P_',71,1;
    'LS_DDI_',73,1;
    'POSTS_DDI_',77,1;
    'LS_AI_',79,1;
    'POSTS_AI_',83,1;
    
    };
for thisArea = 1:size(xlsHeaders,1)
    eval([xlsHeaders{thisArea,1} '= ' num2str(xlsHeaders{thisArea,2}) ';']);
end

%% Define our area types here
global GMTypes;
GMTypes = { % 'Type', ColorCode (RGB);
    'GM',[0.6 0.6 0.6];    % Gray matter without visual modulation
    'VM',[1 0 1];          % Visual modulation but not sure to which area it belongs
    'MST',[0 0 1];
    'MT',[1 0 0];
    'VIP',[1 .6 0];
    'LIP',[0 1 0];
    };
for thisArea = 1:length(GMTypes)
    eval([GMTypes{thisArea,1} '= -' num2str(thisArea) ';']);
end


%% Read Data
[num,txt,raw] = xlsread('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',2);

% Align headers!!!
txt = txt(3:end,:); % txt begins at the 2nd line and ends at the end of data
len = size(txt,1);  % Now it's the real len of data

num = num(4:len+4-1,:);  % num begins at the 1st line
raw = raw(4:len+4-1,:); % raw begins at the 1st line

%% RF plots
% % {

% Pre-allocation
for i = 1:size(xlsHeaders,1)
    if xlsHeaders{i,3} > 0    % Certain length
        eval(['xlsData.' xlsHeaders{i,1} ' = NaN(len,' num2str(xlsHeaders{i,3}) ');']);
    else  % Uncertain length
        eval(['xlsData.' xlsHeaders{i,1} ' = cell(len,1);']);
    end
end


for i = 1 : len
    
    if isnan(raw{i,DEPTH_}); continue; end;  % If there's no depth data, we skip it.
    
    % Decode xls and pack data into xlsData
    xlsData.GRIDLOC_(i,:) = [num(i,LOCX_) num(i,LOCY_)];   % [x,y]
    xlsData.HEMISPHERE_(i,:) = txt{i,4};
    xlsData.DEPTH_(i,:) = raw{i,DEPTH_};      % Get depth
    xlsData.TYPE_(i,:) = GetAreaType(raw{i,TYPE_});      % Get area types
    
    
    if ~isnan(raw{i,RF_}); xlsData.RF_(i,:) = str2num(raw{i,RF_}); end                % Get RF [x,y,w,h]
    
end

% assignin('base','raw',{num,txt,raw});
% assignin('base','xlsData',xlsData);


% Preprocess
uniqueTypes = sort(unique(xlsData.TYPE_(~isnan(xlsData.TYPE_))),'descend');  % AreaTypes
uniqueSites = munique([xlsData.GRIDLOC_ xlsData.DEPTH_]);   % Use gridLoc and depth to define unique recording sites

if ishandle(30); close(30);end;   set(figure(30),'color','w','position',[ 100   100   788   600]);
set(0,'defaultaxesxtickmode','auto');

MT_RFs = [];

for i = 1:sum(uniqueTypes > -100)  % Sure areas
    
    subplot_tight(ceil(sqrt(sum(uniqueTypes > -100))),ceil(sqrt(sum(uniqueTypes > -100))),i,[0.1 0.1]);
    %     subplot_tight(1,3,i,[0.042 0.042]);
    
    set(title(GMTypes{-uniqueTypes(i),1}),'FontSize',15,'Color',GMTypes{-uniqueTypes(i),2});
    axis equal; box on; axis([-45 45 -45 45]);
    line([-45 45],[0 0],'color','k','LineStyle',':'); hold on;
    line([0 0],[-45 45],'color','k','LineStyle',':');
    set(gca,{'xtick','ytick'},{-40:20:40,-40:20:40});
    %     xlabel('degree');
    
    % masks goes here
    sure = find(~isnan(xlsData.RF_(:,1)) & xlsData.TYPE_ == uniqueTypes(i) & (1| num(:, UNITS_)>=5 ))';
    
    % Limit the plotting number for clarity
    %     plotMaxN = 50;
    %     if length(sure) > plotMaxN
    %         sure = sure(randperm(length(sure),plotMaxN));
    %     end
    
    for j = sure   % Sure
        
        % Default: right hemisphere
        if xlsData.HEMISPHERE_(j,:) == 'L'
            xlsData.RF_(j,1) = - xlsData.RF_(j,1);
        end
        
        rectangle('position',[xlsData.RF_(j,1)-xlsData.RF_(j,3)/2 xlsData.RF_(j,2)-xlsData.RF_(j,4)/2, xlsData.RF_(j,3) xlsData.RF_(j,4)],...
            'Curvature',[0.2 0.2],'EdgeColor',GMTypes{-uniqueTypes(i),2},'LineWidth',1.5);
        area(i,j) = xlsData.RF_(j,3)*xlsData.RF_(j,4);
        eccentricity(i,j) = sqrt(xlsData.RF_(j,1)^2+xlsData.RF_(j,2)^2);
        
        if xlsData.TYPE_(j) == MT  % Save MT's location and RF for later plotting
            MT_RFs(end+1,:) = [xlsData.GRIDLOC_(j,:) xlsData.RF_(j,:)];  % [GridX GridY RF_X RF_Y RF_W RF_H]
        end
    end
    
    for j = find(~isnan(xlsData.RF_(:,1)) & xlsData.TYPE_ == uniqueTypes(i) - 100)'   % Not sure
        rectangle('position',[xlsData.RF_(j,1)-xlsData.RF_(j,3)/2 xlsData.RF_(j,2)-xlsData.RF_(j,4)/2, xlsData.RF_(j,3) xlsData.RF_(j,4)],...
            'Curvature',[0.2 0.2],'EdgeColor',GMTypes{-uniqueTypes(i),2},'LineWidth',2,'LineStyle','-');
        area(i,j) = xlsData.RF_(j,3)*xlsData.RF_(j,4);
        eccentricity(i,j) = sqrt(xlsData.RF_(j,1)^2+xlsData.RF_(j,2)^2);
    end
end

% Area Distribution
figure;
c={'b','r','y','g'};

for i = [1 3 2 4]
    hist(sqrt(unique(area(i,area(i,:)~=0))),0:5:100);
    hold on;
    line([mean(sqrt(unique(area(i,area(i,:)~=0)))) mean(sqrt(unique(area(i,area(i,:)~=0))))],[0 15],'linestyle','--','color',c{i},'linewidth',2);
    set(findobj(gca,'type','patch','edgecolor','k'),'edgecolor','none','facecolor',c{i},'facealpha',0.5);
    
end
line([0 80],[.01 .01],'color','k');
line([0.01 .01],[0 15],'color','k');
xlim([0 80]);
xlabel('Equivalent diameter of receptive field (degree)');

% Eccentricity vs Size
figure();
for i = [1 3 2 4]
    plot(eccentricity(i,area(i,:)~=0),sqrt(area(i,area(i,:)~=0)),[c{i} 'o'],'markersize',3,'markerfacecol',c{i}); hold on;
end


% MT: RF v.s. Loc
figure();
plot(MT_RFs(:,1),MT_RFs(:,4),'or','markersize',10,'markerfacecol','r');
xlabel('Grid location X ( A --> P )');
ylabel('MT RF X ( Down --> Up )');
SetFigure(20);
para = polyfit(MT_RFs(:,1),MT_RFs(:,4),1);
hold on; plot([13 25],polyval(para,[13 25]),'r','linewidth',2);



%}

%% Basic Tuning Properties

%%
mask = find(strcmp(txt(:,TYPE_),'MST') & (strcmp(txt(:,Protocol_),'1DAT,NC')|  strcmp(txt(:,Protocol_),'1DAT')) &...
    num(:,UNITS_) == 1 & num(:,TUNING_P_) <= 0.05 );  % 1 MU, 5 SU1, 6 SU2

a = num(mask,PREF_:D270_); % Tuning data
a = munique(a);

[~,ind]=sort(a(:,2),'descend');
a = a(ind,:);

% From heading to polar
a(:,1) = mod(90-a(:,1),360);

set(0,'defaultaxesxtickmode','auto');
figure(11);
clf
roseInHeading(a(:,1)/180*pi);

figure(12);  clf
set(12,'color','w');

% figure; hist(a(:,1),20)
for thisArea = 1: size(a,1)
    polarwitherrorbar([ 0 a(thisArea,1)/180*pi],[0 a(thisArea,2)],[0 0],'b');
    hold on;
end

title(num2str(length(mask)));
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'fontsize',10),'fontsize',14);

%% d90 v.s. d270
figure(13);  clf
hold on
%plot(abs(b(:,5)),abs(b(:,6)),'ro','markersize',10);
plot(abs(a(:,5)),abs(a(:,6)),'bo','markersize',10,'markerfacecolor','b');
axis([0 16 0 16]);
hold on; plot([0 16],[0 16],'k-');

means = mean(abs(a(:,5:6)));
errs = std(abs(a(:,5:6)))/sqrt(size(a,1));

plot(means(1),means(2),'b*','markersize',10);
plot([means(1)-errs(1) means(1)+errs(1)],[means(2) means(2)],'b-','linewidth',2);
plot([means(1) means(1)], [means(2)-errs(2) means(2)+errs(2)],'b-','linewidth',2);

axis square;
[h p ] = ttest(abs(a(:,5))-abs(a(:,6)),0,0.05,'right')
title(['p=' num2str(p)]);
% [h p ] = ttest(abs(b(:,5))-abs(b(:,6)),0,0.05)%,'right')

% d' @ 0 v.s. preferred direction

% figure();
% plot(a(:,1),abs(a(:,5))./abs(a(:,6)),'o'); hold on;
% xlim([0 360]);
% plot([0 360],[0 0],'k-');

%% DDI v.s. HTI

mask = find( (1|strcmp(txt(:,TYPE_),'MST')) & (1|strcmp(txt(:,Protocol_),'1DAT,NC')|  strcmp(txt(:,Protocol_),'1DAT')) &...
    ( 1|num(:,UNITS_) == 1) & (1 | num(:,TUNING_P_) <= 0.05) );  % 1 MU, 5 SU1, 6 SU2

data = num(mask,DDI_:TUNING_P_);

data = munique(data);

figure(14); clf
scatter(data(data(:,3)<0.05,1),data(data(:,3)<0.05,2),7 * -log(data(data(:,3)<0.05,3)+eps),'or');
hold on;
plot(data(data(:,3)>=0.05,1), data(data(:,3)>=0.05,2),'ok');
title(sprintf('N = %g',length(data)));
plot([0 1],[0 1],'--k');
legend('p<0.05','n.s.','location','nw');
axis square;
axis([0 1 0 1]);
xlabel('DDI');
ylabel('HTI');
SetFigure(15);

%% DDI Distribution

clear data;
% areas = {'LIP', 'MST', 'VIP','MT'};
% colors = {'g','b',[1 0.8 0.2],'r'};
areas = {'MST','VIP','MT','LIP' };
colors = {'r','r','r','r'}
figure(20);  clf


for thisArea = 1:length(areas)
    
    mask_S = find(strcmp(txt(:,TYPE_),areas(thisArea)) & (strcmp(txt(:,Protocol_),'1DAT,NC')| strcmp(txt(:,Protocol_),'1DAT')) &...
        (1| num(:,UNITS_) >= 5) & num(:,TUNING_P_) <= 0.05 );  % 1 MU, 5 SU1, 6 SU2
    mask_NS = find(strcmp(txt(:,TYPE_),areas(thisArea)) & (strcmp(txt(:,Protocol_),'1DAT,NC')| strcmp(txt(:,Protocol_),'1DAT')) &...
        (1 | num(:,UNITS_) >= 5) & num(:,TUNING_P_) > 0.05 );  % 1 MU, 5 SU1, 6 SU2
    
    mask_S_SU = find(strcmp(txt(:,TYPE_),areas(thisArea)) & (strcmp(txt(:,Protocol_),'1DAT,NC')| strcmp(txt(:,Protocol_),'1DAT')) &...
        (0| num(:,UNITS_) >= 5) & num(:,TUNING_P_) <= 0.05 );  % 1 MU, 5 SU1, 6 SU2
    mask_NS_SU = find(strcmp(txt(:,TYPE_),areas(thisArea)) & (strcmp(txt(:,Protocol_),'1DAT,NC')| strcmp(txt(:,Protocol_),'1DAT')) &...
        (0| num(:,UNITS_) >= 5) & num(:,TUNING_P_) > 0.05 );  % 1 MU, 5 SU1, 6 SU2
    
    DDI_S = munique(num(mask_S,DDI_));
    DDI_S_SU = munique(num(mask_S_SU,DDI_));
    
    DDI_NS = munique(num(mask_NS,DDI_));
    DDI_NS_SU = munique(num(mask_NS_SU,DDI_));
    
    xcenters = 0:0.05:1;
    
    subplot(2,2,thisArea);
    histS = hist(DDI_S,xcenters);
    histNS = hist(DDI_NS,xcenters);
    hbars = bar(xcenters,[histS' histNS'],1,'stacked','facecolor',colors{thisArea},'LineWidth',2,'edgecolor',colors{thisArea});
    set(hbars(2),'FaceColor','none');
    
    meanDDI = mean([DDI_S; DDI_NS]);
    maxY = max(histS+histNS)*1.4;
    hold on; plot([meanDDI meanDDI],[0 maxY],'--','linewidth',2,'color',colors{thisArea});
    
    SU_S = length(DDI_S_SU);
    SU_NS = length(DDI_NS_SU);
    MU_S = length(DDI_S)-SU_S;
    MU_NS = length(DDI_NS)-SU_NS;
    
    text(0.05, maxY/1.2+0.5, sprintf('Mean = %2.2g\nSU %g/%g = %2.0f%%\nMU %g/%g = %2.0f%%',...
        meanDDI,SU_S, SU_S+SU_NS,SU_S/(SU_S+SU_NS)*100, MU_S, MU_S+MU_NS,MU_S/(MU_S+MU_NS)*100 ));
    
    xlabel('DDI'); ylabel('Number of cases');
    axis([0 1 0 maxY]);
    
    data{thisArea}.DDIall = [DDI_S ; DDI_NS];
    
    SetFigure(15)
    
end;

[~, VIPLIP_DDI_p] = ttest2(data{1}.DDIall,data{3}.DDIall)
VIPLIP_DDI_roc = 1-rocN(data{1}.DDIall,data{3}.DDIall)


%% MU Clustering Analysis

areas = {'VIP','MST'};
stimtypes = {'Vestibular','Vis'};
% colors = {[1 0.8 0.2],'b'};
colors = {'k','k'};

for thisArea = 1:length(areas)
    
    for thisStimtype = 1:length(stimtypes)
        
    MUClustering = nan(1000,3);
    i = 0;
    
        % Overall mask: MU, all tuning
        mask = find(strcmp(txt(:,TYPE_),areas(thisArea)) & (strcmp(txt(:,Protocol_),'1DAT,NC')| strcmp(txt(:,Protocol_),'1DAT')) &...
            ( 0| num(:,UNITS_) == 1) & (1 | num(:,TUNING_P_) <= 0.05) );  % 1 MU, 5 SU1, 6 SU2
        
        session = num(mask,[SESSION_ LOCX_ LOCY_]);
        depth = num(mask,DEPTH_);
        
        if strcmp(stimtypes(thisStimtype),'Vestibular')
            tuningCurves = raw(mask,TUNING_MEAN_VEST_);  % in String
        elseif strcmp(stimtypes(thisStimtype),'Vis')
            tuningCurves = raw(mask,TUNING_MEAN_VIS_);  % in String
        end
        
        [unique_session,counts] = munique(session);
        moreThanOneMUSession = find(counts >= 2)';
        
        for k = moreThanOneMUSession
            % Find lines that belongs to this session
            thisSession = find(sum(abs(session - repmat(unique_session(k,:),size(session,1),1)),2)==0);
            
            thisTuningCurves = zeros(length(thisSession),16);
            
            % Get tuning curves and interpolate to 22.5 interval
            for j = 1:length(thisSession)
                tuningOuter = str2num(tuningCurves{thisSession(j)});
                switch length(tuningOuter)
                    case 8  % 45 degree interval
                        azi = [0:45:359 360];
                        thisTuningCurves(j,:) = interp1(azi,[tuningOuter tuningOuter(1)],0:22.5:359);
                    case 10 % +/- 22.5 added
                        azi = [0 45 67.5 90 112.5 135:45:315 360];
                        thisTuningCurves(j,:) = interp1(azi,[tuningOuter tuningOuter(1)],0:22.5:359);
                    case 16 % 22.5 interval
                        thisTuningCurves(j,:) = tuningOuter;  % do nothing
                end
            end
            
            % For all possible combinations
            for outer = 1:length(thisSession)-1
                depthOuter = depth(thisSession(outer));
                for inner = outer+1 :length(thisSession)
                    depthInner = depth(thisSession(inner));
                    
                    % Tuning correlation
                    [r, p] = corrcoef(thisTuningCurves(outer,:),thisTuningCurves(inner,:));
                    
                    % Save result
                    i = i + 1;
                    MUClustering(i,:) = [abs(depthOuter-depthInner) r(2) p(2)];
                end
            end
        end
        MUClustering(isnan(MUClustering(:,1)),:) = [];
        
%         figure(100+thisArea * length(stimtypes) + thisStimtype); clf ; hist(MUClustering(MUClustering(:,1)<500,1),30); 
%         set(gcf,'name',[areas{thisArea} ', ' stimtypes{thisStimtype}]);
        
        depthDiffGroups = [80 120; 180 220; 280 320; 380 420 ; 420 10000];
        
        figure(200+thisArea * length(stimtypes) + thisStimtype);  clf
        set(gcf,'name',[areas{thisArea} ', ' stimtypes{thisStimtype}]);

        xcenters = -1:0.1:1;
        
        for depthInd = 1:size(depthDiffGroups,1)
            
            r_sig = MUClustering((MUClustering(:,1)>=depthDiffGroups(depthInd,1)) & (MUClustering(:,1)<=depthDiffGroups(depthInd,2))  ...
                & (MUClustering(:,3)<=0.05),2);
            r_ns = MUClustering((MUClustering(:,1)>=depthDiffGroups(depthInd,1)) & (MUClustering(:,1)<=depthDiffGroups(depthInd,2))  ...
                & (MUClustering(:,3)>0.05),2);
            
            histS = hist(r_sig,xcenters);
            histNS = hist(r_ns ,xcenters);
            subplot(size(depthDiffGroups,1),1,depthInd);
            hbars = bar(xcenters,[histS' histNS'],0.8,'stacked','facecolor',colors{thisArea},'LineWidth',2,'edgecolor',colors{thisArea});
            set(hbars(2),'FaceColor','none');
            
            % Mean r test
            
            meanValue = median([r_sig; r_ns]);
            p_meanValue = signtest([r_sig; r_ns],0);
            
            % Annotation ( Portable:) )
            lims = axis(); maxY = lims(4);
            arrowSize = 0.08 * [1 5*abs((lims(4)-lims(3))/(lims(2)-lims(1)))]; % Keep it looking good
            hold on;
            fill([meanValue meanValue+arrowSize(1)/2 meanValue-arrowSize(1)/2],[maxY+arrowSize(2) maxY+arrowSize(2)*2 maxY+arrowSize(2)*2],'k');
            plot([meanValue meanValue],[maxY+arrowSize(2)*2 maxY+arrowSize(2)*3],'k');
            ylim([lims(3) maxY+arrowSize(2)*3]);
            text(0, maxY+arrowSize(2)*3, sprintf('Median = %2.2f\n\\itp \\rm= %.2g',meanValue,p_meanValue),'HorizontalAlignment','right');
            
            xlim([-1.05 1.05]);
            
            if depthInd == ceil(size(depthDiffGroups,1)/2)
                ylabel('Number of cases');
            end
            
            if depthInd ~= size(depthDiffGroups,1)
                set(gca,'xticklabe',[]);
            else
                xlabel('Correlation coefficient');
            end
            
        end
        
        SetFigure(15);
        
    end;
end;

%% LIP and VIP: Tuning and Delayed Saccade Analysis
clear data;

areas = {'LIP','VIP'};
colors = {'g',[1 0.8 0.2]};
% colors = {'k','k'};
figure(22);  clf

for thisArea = 1:length(areas)
    
    % Overall mask: MU, all tuning
    mask = find(strcmp(txt(:,TYPE_),areas(thisArea)) & strcmp(txt(:,Protocol_),'DelSac') &...
        ( 1 | num(:,UNITS_) == 1) & (1 | num(:,TUNING_P_) <= 0.05) );  % 1 MU, 5 SU1, 6 SU2
    
    data{thisArea}.dsP = num(mask,LS_P_:POSTS_P_);
    data{thisArea}.dsP(:,6) = num(mask,TUNING_P_);
    data{thisArea}.dsDDI = num(mask,LS_DDI_:POSTS_DDI_);
    data{thisArea}.dsDDI(:,6) = num(mask,DDI_);
    data{thisArea}.dsAI = num(mask,LS_AI_:POSTS_AI_);
    
    xcenters = 0:0.05:1;
    
    for phase = 1:5  %  LS Sustain Pre Co Post 1-DAzimuthTuning
        
        DDI_sig = data{thisArea}.dsDDI(data{thisArea}.dsP(:,phase)<=0.05,phase);
        DDI_ns = data{thisArea}.dsDDI(data{thisArea}.dsP(:,phase) > 0.05,phase);
        
        histS = hist(DDI_sig,xcenters);
        histNS = hist(DDI_ns ,xcenters);
        
        subplot(2,5,(thisArea -1) * 5 + phase);
        
        hbars = bar(xcenters,[histS' histNS'],1,'stacked','facecolor',colors{thisArea},'LineWidth',2,'edgecolor',colors{thisArea});
        set(hbars(2),'FaceColor','none');
        
        % Mean r test
        
        meanValue = median([DDI_sig; DDI_ns]);
        % p_meanValue = signtest([DDI_sig; DDI_ns],0);
        
        % Annotation ( Portable:) )
        lims = axis(); maxY = lims(4);
        arrowSize = 0.05 * [1 5*abs((lims(4)-lims(3))/(lims(2)-lims(1)))]; % Keep it looking good
        hold on;
        fill([meanValue meanValue+arrowSize(1)/2 meanValue-arrowSize(1)/2],[maxY+arrowSize(2) maxY+arrowSize(2)*2 maxY+arrowSize(2)*2],'k');
        plot([meanValue meanValue],[maxY+arrowSize(2)*2 maxY+arrowSize(2)*3],'k');
        ylim([lims(3) maxY+arrowSize(2)*3]);
        text(meanValue+arrowSize(1)*2, maxY+arrowSize(2)*3, sprintf('Mean = %2.2f\n%g/%g = %2.2g%%',...
            meanValue,length(DDI_sig),length(DDI_sig)+length(DDI_ns),length(DDI_sig)/(length(DDI_sig)+length(DDI_ns))*100));
        
        xlim([min(xcenters) max(xcenters)]);
        
        
    end
    
end

% ttest and ROC analysis pooled data

for phase = 1:5  %  LS Sustain Pre Co Post 1-DAzimuthTuning
    [~,pp(phase)] = ttest2(data{1}.dsDDI(:,phase),data{2}.dsDDI(:,phase));
    roc(phase) = rocN(data{1}.dsDDI(:,phase),data{2}.dsDDI(:,phase));
end

disp(pp)
disp(roc)

%% PCA analysis

% Pack data matrix
index = 1;
dataPCA = [];

for thisArea = 1:length(areas)
    dataPCA = [dataPCA; thisArea *ones(size(data{thisArea}.dsDDI,1),1) data{thisArea}.dsP <= 0.05 data{thisArea}.dsDDI];
end
dataPCA(isnan([data{1}.dsDDI(:,end); data{2}.dsDDI(:,end)]),:) = [];

% Do PCA
[sortedEigVectors, ~] = pca(dataPCA(:,2:7));

% The first three eigenvectors that have the first three largest
% eigenvalues.
PC1 = sortedEigVectors(:,1);
PC2 = sortedEigVectors(:,2);
% PC3 = sortedEigVectors(:,3);

% Projecting the raw data onto the first three eigenvectors
projPC1 = dataPCA(:,2:7) * PC1;
projPC2 = dataPCA(:,2:7) * PC2;
% projPC3 = dataPCA(:,2:end) * PC3;

figure(23); clf
for thisArea = 1:length(areas)
    %     plot3(projPC1(dataPCA(:,1)==thisArea), projPC2(dataPCA(:,1)==thisArea), projPC3(dataPCA(:,1)==thisArea),...
    %         'linestyle','none','marker','o','markerFaceColor',colors{thisArea},'markerEdgeColor',colors{thisArea}); grid on; hold on;
    %     xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    
    plot(projPC1(dataPCA(:,1)==thisArea), projPC2(dataPCA(:,1)==thisArea), ...
        'linestyle','none','marker','o','markerFaceColor',colors{thisArea},'markerEdgeColor',colors{thisArea}); grid on; hold on;
    xlabel('PC1'); ylabel('PC2');
    
end

figure(24); clf
subplot(3,1,1); bar(PC1);
subplot(3,1,2); bar(PC2);
% subplot(3,1,3); bar(PC3);

%% Direct 2-D plot
figure(25); clf
xx = dataPCA(:,13); % Tuning DDI;
yy = mean(dataPCA(:,12),2); % average DS DDI (2-5);

for thisArea = 1:length(areas)
    plot(xx(dataPCA(:,1)==thisArea), yy(dataPCA(:,1)==thisArea), ...
        'linestyle','none','marker','o','markerFaceColor',colors{thisArea},'markerEdgeColor',colors{thisArea}); grid on; hold on;
    xlabel('1-D Tuning DDI'); ylabel('Del Sac DDI');
end



