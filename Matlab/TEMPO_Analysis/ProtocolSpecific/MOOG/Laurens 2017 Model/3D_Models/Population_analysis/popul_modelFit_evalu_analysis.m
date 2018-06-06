% population analysis of those 3D models
% VAJ, VA, VJ, AJ, VO, AO
% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
% LBY 20170406
% LBY 20170605

% function popul_modelFit_evalu_analysis(models,Protocol)

%% load data
clear all;
models = {'VO','AO','VA','VAJ'};
% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
Protocol = 1;
FnameCode = 2;



switch Protocol
    case 1
        cd('Z:\Data\TEMPO\BATCH\PCC_QQ_PSTH3DModel');
        % load('PSTH3DModel_T_OriData.mat');
        load('PSTH3DModel_T_PARAs.mat');
    case 2
        cd('Z:\Data\TEMPO\BATCH\PCC_QQ_PSTH3DModel_R');
        % load('PSTH3DModel_R_OriData.mat');
        load('PSTH3DModel_R_PARAs.mat');
end
disp('Load PARAs data SUCCESS!');

%% analysis for 3D model parameters
% BIC, RSS, R_squared and weight

unique_elevation = [-90 -45 0 45 90];
unique_azimuth = [0 45 90 135 180 225 270 315];

% Compare the BIC, RSS and R_squared value across models
% to find the best

for cell_inx = 1:length(CellName)
    for k = 1:length(StimType{cell_inx})
        for m_inx = 1:length(models)
            
            eval(['BIC_Norm{', num2str(k), '}(', num2str(cell_inx), ',',num2str(m_inx), ')= BIC(', num2str(cell_inx), ').', models{m_inx}, '{', num2str(k), '};']);
            eval(['RSS_Norm{', num2str(k), '}(', num2str(cell_inx), ',',num2str(m_inx),  ')= RSS(', num2str(cell_inx), ').', models{m_inx}, '{', num2str(k), '};']);
            eval(['R_squared_Norm{', num2str(k), '}(', num2str(cell_inx), ',',num2str(m_inx),  ')= R_squared(', num2str(cell_inx), ').', models{m_inx}, '{', num2str(k), '};']);
            
        end
        
        % find the best fit model, return the index
        [~,BIC_min_inx(cell_inx,k)] = min(BIC_Norm{k}(cell_inx,:));
        [~,RSS_min_inx(cell_inx,k)] = min(RSS_Norm{k}(cell_inx,:));
        [~,R_squared_max_inx(cell_inx,k)] = max(R_squared_Norm{k}(cell_inx,:));
        
    end
end

% %% weights analysis & figures
% 
% % weight for VAJ model
% % w_v, w_a, w_j
% 
% for cell_inx = 1:length(CellName)
%     for k = 1:length(StimType{cell_inx})
%         weight_VAJ{k}(cell_inx,1) = PARA(cell_inx).VAJ{1,k}(1,29)*(1-PARA(cell_inx).VAJ{1,k}(1,30));
%         weight_VAJ{k}(cell_inx,2) = (1-PARA(cell_inx).VAJ{1,k}(1,29))*(1-PARA(cell_inx).VAJ{1,k}(1,30));
%         weight_VAJ{k}(cell_inx,3) = PARA(cell_inx).VAJ{1,k}(1,30);
%     end
% end
% 
% pro = {'Tranlation','Rotation'};
% stimtype = {'vestibular','visual','combined'};
% figure(131);
% set(gcf,'pos',[0 0 1900 1000]);
% clf;
% 
% for k = 1:length(stimtype)
%     subplot(1,3,k);
%     bar(squeeze(weight_VAJ{k}(:,:)),'stack');
%     legend('w v','w a', 'w j');
%     title(['Weights distributions of VAJ model  ',pro(Protocol), '   ',stimtype{k}]);
%     ss = ['Weights distributions of VAJ model',stimtype{k}];
% end
% saveas(131,['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning_models\Population_results\',ss],'emf');
% 
% % weight for VA model
% % w_v, w_a
% 
% for cell_inx = 1:length(CellName)
%     for k = 1:length(StimType{cell_inx})
%         weight_VA{k}(cell_inx,1) = PARA(cell_inx).VA{1,k}(1,21);
%         weight_VA{k}(cell_inx,2) = 1-PARA(cell_inx).VA{1,k}(1,21);
%     end
% end
% 
% pro = {'Tranlation','Rotation'};
% stimtype = {'vestibular','visual','combined'};
% figure(141);
% set(gcf,'pos',[0 0 1900 1000]);
% clf;
% 
% for k = 1:length(stimtype)
%     subplot(1,3,k);
%     bar(squeeze(weight_VA{k}(:,:)),'stack');
%     legend('w v','w a');
%     title(['Weights distributions of VA model  ',pro(Protocol), '   ',stimtype{k}]);
%     ss = ['Weights distributions of VA model',stimtype{k}];
% end
% set(gcf,'paperpositionmode','auto');
% saveas(141,['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning_models\Population_results\',ss],'emf');
%% plot BIC values distributions for each cell of each stimtype

pro = {'Tranlation','Rotation'};
stimtype = {'vestibular','visual','combined'};

for cell_inx = 1:length(CellName)
    figure(151+cell_inx);
    clf;
    set(gcf,'pos',[0 0 1800 500]);
    for k = 1:length(StimType{cell_inx})
        subplot(1,3,k);
        bar(BIC_Norm{k}(cell_inx,:));
        set(gca, 'xticklabel',{'VO','AO','VA','VAJ'});
        title(stimtype{k});
    end
    suptitle(['BIC distributions of ',CellName{cell_inx}]);
    ss = ['BIC_distributions_',CellName{cell_inx}];
    
    set(gcf,'paperpositionmode','auto');
    switch Protocol
        case 1
            saveas(151+cell_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning_models\Population_results\',ss],'emf');
        case 2
            saveas(151+cell_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_R_Tuning_models\Population_results\',ss],'emf');
    end
end


%% plot R_squared distributions & BIC_RSS min value distributions

pro = {'Tranlation','Rotation'};
stimtype = {'vestibular','visual','combined'};

for k = 1:length(stimtype)
    figure(111+k);
    set(gcf,'pos',[0 0 1900 1000]);
    clf;
    
    for m_inx = 1:length(models)
        subplot(2,4,m_inx);
        hist(R_squared_Norm{k}(:,m_inx));
        title(models{m_inx});
        xlim([0 1]);
    end
    
    ss = ['R_squared_distributions_',stimtype{k}];
    suptitle(['R squared distributions     ',pro(Protocol),'   ', stimtype{k}]);
    
    % save figures
    set(gcf,'paperpositionmode','auto');
    switch Protocol
        case 1
            saveas(111+k,['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning_models\Population_results\',ss],'emf');
        case 2
            saveas(111+k,['Z:\LBY\Recording data\Qiaoqiao\3D_R_Tuning_models\Population_results\',ss],'emf');
    end
end

figure(121);
set(gcf,'pos',[0 0 1900 1000]);
clf;
for k = 1:length(stimtype)
    subplot(3,3,1+(k-1)*3);
    hist(BIC_min_inx(:,k));
    xlim([1 7]);
    title('BIC (min) ');
    subplot(3,3,2+(k-1)*3);
    hist(RSS_min_inx(:,k));
    xlim([1 7]);
    title('RSS (min) ');
    subplot(3,3,3+(k-1)*3);
    hist(R_squared_max_inx(:,k));
    xlim([1 7]);
    title('R squared (max)');
    suptitle(['Best models fitting distributions  ',pro(Protocol),'   ', stimtype{k}]);
    ss = ['Best models fitting distributions',stimtype{k}];
    
    % save figures
    set(gcf,'paperpositionmode','auto');
    switch Protocol
        case 1
            saveas(121,['Z:\LBY\Recording data\Qiaoqiao\3D_T_Tuning_models\Population_results\',ss],'emf');
        case 2
            saveas(121,['Z:\LBY\Recording data\Qiaoqiao\3D_R_Tuning_models\Population_results\',ss],'emf');
    end
end

% save the data
switch Protocol
    case 1
        save('PSTH3DBestFitModel_T_PARAs.mat', 'CellName', 'PARA','BIC_Norm','BIC_min_inx','RSS_Norm','RSS_min_inx','R_squared_Norm','R_squared_max_inx','weight_VAJ','weight_VA');
    case 2
        save('PSTH3DBestFitModel_R_PARAs.mat', 'CellName', 'PARA','BIC_Norm','BIC_min_inx','RSS_Norm','RSS_min_inx','R_squared_Norm','R_squared_max_inx','weight_VAJ','weight_VA');
end
disp('DATA SAVED!');

% end


