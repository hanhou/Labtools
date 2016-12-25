
%% Definition. Temporality hard coded
global xlsHeaders;
xlsHeaders = {  % Name, Loc in xls, length
                'GRIDLOC_', 3, 2;
                'DEPTH_', 5, 1;
                'TYPE_', 6, 1;
                'SNR_', 12, 1;
                'SPKWAVE_', 13, -1;
                'RF_', 14, 4;
                'PREF_', 15, 1;
                'DDI_', 16, 1;
                'HTI_', 17, 1;
                'TUNING_P_', 18,1;
                'D90_',19,1;
                'D270_', 20, 1;
                'TUNING_MEAN_', 21, -1;
                'TUNING_ERR_', 22, -1;
        };
for i = 1:size(xlsHeaders,1)
    eval([xlsHeaders{i,1} '= ' num2str(xlsHeaders{i,2}) ';']);  
end

%% Define our area types here
global GMTypes;
GMTypes = { % 'Type', ColorCode (RGB);
    'GM',[0.6 0.6 0.6];    % Gray matter without visual modulation
    'VM',[1 0 1];          % Visual modulation but not sure to which area it belongs
    'MST',[0 0 1]; 
    'MT',[1 0 0];
    'LIP',[0 1 0];
    'VIP',[1 .6 0];
    };
for i = 1:length(GMTypes)
    eval([GMTypes{i,1} '= -' num2str(i) ';']); 
end 


%% Read Data
[num,txt,raw] = xlsread('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',2);
len = size(raw,1);

% Pre-allocation
for i = 1:size(xlsHeaders,1)
    if xlsHeaders{i,3} > 0    % Certain length
        eval(['xlsData.' xlsHeaders{i,1} ' = NaN(len,' num2str(xlsHeaders{i,3}) ');']);
    else  % Uncertain length
        eval(['xlsData.' xlsHeaders{i,1} ' = cell(len,1);']);
        
    end
end

% Process 
for i = 3 : len
   
    % Decode xls and pack data into xlsData
    xlsData.GRIDLOC_(i,:) = str2num(raw{i,GRIDLOC_}(2:end));   % [x,y]
%     xlsData.DEPTH_(i,:) = round(raw{i,DEPTH_}/100);      % Get depth (precision = 100 um)
    xlsData.TYPE_(i,:) = GetAreaType(raw{i,TYPE_});      % Get area types
    if ~isnan(raw{i,RF_}); xlsData.RF_(i,:) = str2num(raw{i,RF_}); end                % Get RF [x,y,w,h]
    
    
end

assignin('base','raw',{num,txt,raw});
assignin('base','xlsData',xlsData);

%% RF plots

% Preprocess
uniqueTypes = sort(unique(xlsData.TYPE_(~isnan(xlsData.TYPE_))),'descend');  % AreaTypes
uniqueSites = munique([xlsData.GRIDLOC_ xlsData.DEPTH_]);   % Use gridLoc and depth to define unique recording sites

if ishandle(30); close(30);end;   set(figure(30),'color','w','position',[ 560   184   788   764]);
set(0,'defaultaxesxtickmode','auto');

MT_RFs = [];

for i = 1:sum(uniqueTypes > -100)  % Sure areas
    
    %subplot_tight(ceil(sqrt(sum(uniqueTypes > -100))),ceil(sqrt(sum(uniqueTypes > -100))),i,[0.07 0.07]);
    subplot_tight(1,3,i,[0.042 0.042]);
    
    set(title(GMTypes{-uniqueTypes(i),1}),'FontSize',15,'Color',GMTypes{-uniqueTypes(i),2});
    axis equal; box on; axis([-45 45 -45 45]);
    line([-45 45],[0 0],'color','k','LineStyle',':'); hold on;
    line([0 0],[-45 45],'color','k','LineStyle',':');
    set(gca,{'xtick','ytick'},{-40:20:40,-40:20:40});
    xlabel('degree');

    for j = find(~isnan(xlsData.RF_(:,1)) & xlsData.TYPE_ == uniqueTypes(i))'   % Sure
        rectangle('position',[xlsData.RF_(j,1)-xlsData.RF_(j,3)/2 xlsData.RF_(j,2)-xlsData.RF_(j,4)/2, xlsData.RF_(j,3) xlsData.RF_(j,4)],...
                    'Curvature',[0.2 0.2],'EdgeColor',GMTypes{-uniqueTypes(i),2},'LineWidth',1.5);
        area(i,j) = xlsData.RF_(j,3)*xlsData.RF_(j,4);
        
        if xlsData.TYPE_(j) == MT  % Save MT's location and RF for later plotting
            MT_RFs(end+1,:) = [xlsData.GRIDLOC_(j,:) xlsData.RF_(j,:)];  % [GridX GridY RF_X RF_Y RF_W RF_H]
        end
    end
   
    for j = find(~isnan(xlsData.RF_(:,1)) & xlsData.TYPE_ == uniqueTypes(i) - 100)'   % Not sure
        rectangle('position',[xlsData.RF_(j,1)-xlsData.RF_(j,3)/2 xlsData.RF_(j,2)-xlsData.RF_(j,4)/2, xlsData.RF_(j,3) xlsData.RF_(j,4)],...
                    'Curvature',[0.2 0.2],'EdgeColor',GMTypes{-uniqueTypes(i),2},'LineWidth',2,'LineStyle','-');
        area(i,j) = xlsData.RF_(j,3)*xlsData.RF_(j,4);
    end
end

%% Area Distribution
figure;
c={'b','r','g'};

for i = [1 3 2]
    hist(sqrt(area(i,area(i,:)~=0)),0:5:100); hold on;
    line([mean(sqrt(area(i,area(i,:)~=0))) mean(sqrt(area(i,area(i,:)~=0)))],[0 15],'linestyle','--','color',c{i},'linewidth',2);
    set(findobj(gca,'type','patch','edgecolor','k'),'edgecolor','none','facecolor',c{i},'facealpha',0.5);
    
end
line([0 80],[.01 .01],'color','k');
line([0.01 .01],[0 15],'color','k');
xlim([0 80]);
xlabel('Equivalent diameter of receptive field (degree)');


%% MT: RF v.s. Loc
figure();
plot(MT_RFs(:,1),MT_RFs(:,4),'or','markersize',10,'markerfacecol','r');
xlabel('Grid location X ( A --> P )');
ylabel('MT RF X ( Down --> Up )');
SetFigure(20);
para = polyfit(MT_RFs(:,1),MT_RFs(:,4),1);
hold on; plot([13 25],polyval(para,[13 25]),'r','linewidth',2);



