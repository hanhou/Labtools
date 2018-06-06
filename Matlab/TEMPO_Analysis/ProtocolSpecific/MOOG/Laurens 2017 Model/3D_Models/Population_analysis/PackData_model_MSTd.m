% population analysis of those 3D models
% VAJ, VA, VJ, AJ, VO, AO
% Protocol: 1->translation, 2-> rotation
% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
% LBY 20171130


%% load data & pack data
clear all;
cd('Z:\Data\TEMPO\BATCH\MSTd_3DTuning');
load('Z:\Data\TEMPO\BATCH\MSTd_3DTuning\PSTH_OriData.mat');
Monkey = 'Que';

% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
% models = {'VO','AO','VA','VJ','AJ','VAJ'};
% models = {'VO','AO','VA','VAJ'};
models = {'VA'};
global PSTH3Dmodel PSTH;

%%%%%%%%%%%%%%%%%%%%%%%%%% for Translation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_vestiNo = 0;
T_visNo = 0;
for cell_inx = 1:length(QQ_3DTuning_T)
    T_model(cell_inx).name = QQ_3DTuning_T(cell_inx).name;
    T_model(cell_inx).ch = QQ_3DTuning_T(cell_inx).ch;
    
    % for Translation - Vestibular
    if ~isempty(find(QQ_3DTuning_T(cell_inx).stimType == 1))
        T_model(cell_inx).vestiPSTH3Dmodel = QQ_3DTuning_T(cell_inx).PSTH3Dmodel{find(QQ_3DTuning_T(cell_inx).stimType == 1)};
    else
        T_model(cell_inx).vestiPSTH3Dmodel = nan;
    end
    % for Translation - Visual
    if ~isempty(find(QQ_3DTuning_T(cell_inx).stimType == 2))
        T_model(cell_inx).visPSTH3Dmodel = QQ_3DTuning_T(cell_inx).PSTH3Dmodel{find(QQ_3DTuning_T(cell_inx).stimType == 2)};
    else
        T_model(cell_inx).visPSTH3Dmodel = nan;
    end
    
    if isstruct(T_model(cell_inx).vestiPSTH3Dmodel)
        T_vestiNo = T_vestiNo+1;
        T_model(cell_inx).vestiSig = 1;
        % Partial R_squared
        T_PartR2_VAJ_vesti(cell_inx).R2V = T_model(cell_inx).vestiPSTH3Dmodel.VAJ_R2V;
        T_PartR2_VAJ_vesti(cell_inx).R2A = T_model(cell_inx).vestiPSTH3Dmodel.VAJ_R2A;
        T_PartR2_VAJ_vesti(cell_inx).R2J = T_model(cell_inx).vestiPSTH3Dmodel.VAJ_R2J;
        T_wVAJ_vesti(cell_inx).wV =  T_model(cell_inx).vestiPSTH3Dmodel.VAJ_wV;
        T_wVAJ_vesti(cell_inx).wA =  T_model(cell_inx).vestiPSTH3Dmodel.VAJ_wA;
        T_wVAJ_vesti(cell_inx).wJ =  T_model(cell_inx).vestiPSTH3Dmodel.VAJ_wJ;
        T_wVA_vesti(cell_inx).wV =  T_model(cell_inx).vestiPSTH3Dmodel.VA_wV;
        T_wVA_vesti(cell_inx).wA =  T_model(cell_inx).vestiPSTH3Dmodel.VA_wA;
        T_PartR2VA_vesti(cell_inx).wV =  T_model(cell_inx).vestiPSTH3Dmodel.VA_R2V;
        T_PartR2VA_vesti(cell_inx).wA =  T_model(cell_inx).vestiPSTH3Dmodel.VA_R2A;
        
        for m_inx = 1:length(models)
            % pack RSS values to RSS.*(* the model)
            eval(['T_RSS_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.rss_', models{m_inx} , ';']);
            
            %             pack R_squared values to RSS.*(* the model)
            eval(['T_Rsquared_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.RSquared_', models{m_inx} , ';']);
            
            % pack BIC values to RSS.*(* the model)
            eval(['T_BIC_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.BIC_', models{m_inx}, ';']);
            
            % pack model fitting parameters to PARA.*.#
            eval(['T_PARA_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.modelFitPara_', models{m_inx}, ';']);
            
        end
    else
        T_model(cell_inx).vestiSig = 0;
        % Partial R_squared
        T_PartR2_VAJ_vesti(cell_inx).R2V = nan;
        T_PartR2_VAJ_vesti(cell_inx).R2A = nan;
        T_PartR2_VAJ_vesti(cell_inx).R2J = nan;
        T_wVAJ_vesti(cell_inx).wV = nan;
        T_wVAJ_vesti(cell_inx).wA = nan;
        T_wVAJ_vesti(cell_inx).wJ = nan;
        T_wVA_vesti(cell_inx).wV = nan;
        T_wVA_vesti(cell_inx).wA = nan;
        T_PartR2VA_vesti(cell_inx).wV =  nan;
        T_PartR2VA_vesti(cell_inx).wA =  nan;
        
        for m_inx = 1:length(models)
            % pack RSS values to RSS.*(* the model)
            eval(['T_RSS_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            % pack R_squared values to RSS.*(* the model)
            eval(['T_Rsquared_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            % pack BIC values to RSS.*(* the model)
            eval(['T_BIC_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            % pack model fitting parameters to PARA.*.#
            eval(['T_PARA_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            
        end
    end
    if isstruct(T_model(cell_inx).visPSTH3Dmodel)
        T_model(cell_inx).visSig = 1;
        T_visNo = T_visNo+1;
        % Partial R_squared
        T_PartR2_VAJ_vis(cell_inx).R2V = T_model(cell_inx).visPSTH3Dmodel.VAJ_R2V;
        T_PartR2_VAJ_vis(cell_inx).R2A = T_model(cell_inx).visPSTH3Dmodel.VAJ_R2A;
        T_PartR2_VAJ_vis(cell_inx).R2J = T_model(cell_inx).visPSTH3Dmodel.VAJ_R2J;
        T_wVAJ_vis(cell_inx).wV =  T_model(cell_inx).visPSTH3Dmodel.VAJ_wV;
        T_wVAJ_vis(cell_inx).wA =  T_model(cell_inx).visPSTH3Dmodel.VAJ_wA;
        T_wVAJ_vis(cell_inx).wJ =  T_model(cell_inx).visPSTH3Dmodel.VAJ_wJ;
        T_wVA_vis(cell_inx).wV =  T_model(cell_inx).visPSTH3Dmodel.VA_wV;
        T_wVA_vis(cell_inx).wA =  T_model(cell_inx).visPSTH3Dmodel.VA_wA;
        T_PartR2VA_vis(cell_inx).wV =  T_model(cell_inx).visPSTH3Dmodel.VA_R2V;
        T_PartR2VA_vis(cell_inx).wA =  T_model(cell_inx).visPSTH3Dmodel.VA_R2A;
        
        for m_inx = 1:length(models)
            % pack RSS values to RSS.*(* the model)
            eval(['T_RSS_vis(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').visPSTH3Dmodel.rss_', models{m_inx} , ';']);
            
            % pack R_squared values to RSS.*(* the model)
            eval(['T_Rsquared_vis(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').visPSTH3Dmodel.RSquared_', models{m_inx} , ';']);
            
            % pack BIC values to RSS.*(* the model)
            eval(['T_BIC_vis(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').visPSTH3Dmodel.BIC_', models{m_inx}, ';']);
            
            % pack model fitting parameters to PARA.*.#
            eval(['T_PARA_vis(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').visPSTH3Dmodel.modelFitPara_', models{m_inx}, ';']);
            
        end
    else
        T_model(cell_inx).visSig = 0;
        % Partial R_squared
        T_PartR2_VAJ_vis(cell_inx).R2V = nan;
        T_PartR2_VAJ_vis(cell_inx).R2A = nan;
        T_PartR2_VAJ_vis(cell_inx).R2J = nan;
        T_wVAJ_vis(cell_inx).wV = nan;
        T_wVAJ_vis(cell_inx).wA = nan;
        T_wVAJ_vis(cell_inx).wJ = nan;
        T_wVA_vis(cell_inx).wV = nan;
        T_wVA_vis(cell_inx).wA = nan;
        T_PartR2VA_vis(cell_inx).wV =  nan;
        T_PartR2VA_vis(cell_inx).wA =  nan;
        
        for m_inx = 1:length(models)
            % pack RSS values to RSS.*(* the model)
            eval(['T_RSS_vis(',num2str(cell_inx),').', models{m_inx}, ' = nan;']);
            
            % pack R_squared values to RSS.*(* the model)
            eval(['T_Rsquared_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            % pack BIC values to RSS.*(* the model)
            eval(['T_BIC_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
            % pack model fitting parameters to PARA.*.#
            eval(['T_PARA_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
            
        end
    end
end

disp('T model data loaded SUCCESS!');



% save the data
save('PSTH3DModel_T_OriData.mat','T_model','T_PARA_vis', 'T_RSS_vis', 'T_BIC_vis','T_Rsquared_vis','T_PARA_vesti', 'T_RSS_vesti', 'T_BIC_vesti','T_Rsquared_vesti','T_vestiNo','T_visNo','T_PartR2_VAJ_vis','T_PartR2_VAJ_vesti','T_PartR2VA_vis','T_PartR2VA_vis','T_wVAJ_vis','T_wVAJ_vesti','T_wVA_vis','T_wVA_vesti');
disp('DATA SAVED!');


%             w_paras = {{'wV'},...
%                 {'wA'},...
%                 {'wV','wA'},...
%                 {'wV','wJ'},...
%                 {'wA','wJ'},...
%                 {'wV','wA','wJ'},...
%                 {'wV','wA','wJ','wP'},...
%                 };
