% population analysis of those 3D models
% VAJ, VA, VJ, AJ, VO, AO
% Protocol: 1->translation, 2-> rotation
% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
% LBY 20171130


%% load data & pack data
clear all;
cd('Z:\Data\TEMPO\BATCH\FEFp_3DTuning');
load('Z:\Data\TEMPO\BATCH\FEFp_3DTuning\PSTH_OriData.mat');
Monkey = 'Que';

% models = {'VO','AO','VA','VJ','AJ','VAJ','PVAJ'};
models = {'VO','AO','VA','VJ','AJ','VAJ'};
% models = {'VO','AO','VA','VAJ'};
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
    
    
    if isstruct(T_model(cell_inx).vestiPSTH3Dmodel)
        T_vestiNo = T_vestiNo+1;
        T_model(cell_inx).vestiSig = 1;
        % Partial R_squared
        T_PartR2_vesti(cell_inx).R2V = T_model(cell_inx).vestiPSTH3Dmodel.R2V;
        T_PartR2_vesti(cell_inx).R2A = T_model(cell_inx).vestiPSTH3Dmodel.R2A;
        T_PartR2_vesti(cell_inx).R2J = T_model(cell_inx).vestiPSTH3Dmodel.R2J;
        T_wVAJ_vesti(cell_inx).wV =  T_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(29)*(1-T_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30));
        T_wVAJ_vesti(cell_inx).wA =  (1-T_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(29))*(1-T_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30));
        T_wVAJ_vesti(cell_inx).wJ =  T_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30);
        
        for m_inx = 1:length(models)
            % pack RSS values to RSS.*(* the model)
            eval(['T_RSS_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.rss_', models{m_inx} , ';']);
            
            % pack R_squared values to RSS.*(* the model)
            eval(['T_Rsquared_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.RSquared_', models{m_inx} , ';']);
            
            % pack BIC values to RSS.*(* the model)
            eval(['T_BIC_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.BIC_', models{m_inx}, ';']);
            
            % pack model fitting parameters to PARA.*.#
            eval(['T_PARA_vesti(',num2str(cell_inx),').', models{m_inx},' = T_model(',num2str(cell_inx),').vestiPSTH3Dmodel.modelFitPara_', models{m_inx}, ';']);
            
        end
    else
        T_model(cell_inx).vestiSig = 0;
        % Partial R_squared
        T_PartR2_vesti(cell_inx).R2V = nan;
        T_PartR2_vesti(cell_inx).R2A = nan;
        T_PartR2_vesti(cell_inx).R2J = nan;
        T_wVAJ_vesti(cell_inx).wV = nan;
        T_wVAJ_vesti(cell_inx).wA = nan;
        T_wVAJ_vesti(cell_inx).wJ = nan;
        
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
    
end

disp('T model data loaded SUCCESS!');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% for Rotation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_vestiNo = 0;
% R_visNo = 0;
% for cell_inx = 1:length(QQ_3DTuning_R)
%     R_model(cell_inx).name = QQ_3DTuning_R(cell_inx).name;
%     R_model(cell_inx).ch = QQ_3DTuning_R(cell_inx).ch;
%     % for Rotation - Vestibular
%     
%     if ~isempty(find(QQ_3DTuning_R(cell_inx).stimType == 1))
%         R_model(cell_inx).vestiPSTH3Dmodel = QQ_3DTuning_R(cell_inx).PSTH3Dmodel{find(QQ_3DTuning_R(cell_inx).stimType == 1)};
%     else
%         R_model(cell_inx).vestiPSTH3Dmodel = nan;
%     end
%     
%     % for Rotation - Visual
%     if ~isempty(find(QQ_3DTuning_R(cell_inx).stimType == 2))
%         R_model(cell_inx).visPSTH3Dmodel = QQ_3DTuning_R(cell_inx).PSTH3Dmodel{find(QQ_3DTuning_R(cell_inx).stimType == 2)};
%     else
%         R_model(cell_inx).visPSTH3Dmodel = nan;
%     end
%     
%     if isstruct(R_model(cell_inx).vestiPSTH3Dmodel)
%         R_model(cell_inx).vestiSig = 1;
%         R_vestiNo = R_vestiNo+1;
%         % Partial R_squared
%         R_PartR2_vesti(cell_inx).R2V = R_model(cell_inx).vestiPSTH3Dmodel.R2V;
%         R_PartR2_vesti(cell_inx).R2A = R_model(cell_inx).vestiPSTH3Dmodel.R2A;
%         R_PartR2_vesti(cell_inx).R2J = R_model(cell_inx).vestiPSTH3Dmodel.R2J;
%         R_wVAJ_vesti(cell_inx).wV =  R_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(29)*(1-R_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30));
%         R_wVAJ_vesti(cell_inx).wA =  (1-R_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(29))*(1-R_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30));
%         R_wVAJ_vesti(cell_inx).wJ =  R_model(cell_inx).vestiPSTH3Dmodel.modelFitPara_VAJ(30);
%         
%         for m_inx = 1:length(models)
%             % pack RSS values to RSS.*(* the model)
%             eval(['R_RSS_vesti(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').vestiPSTH3Dmodel.rss_', models{m_inx} , ';']);
%             
%             % pack R_squared values to RSS.*(* the model)
%             eval(['R_Rsquared_vesti(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').vestiPSTH3Dmodel.RSquared_', models{m_inx} , ';']);
%             
%             % pack BIC values to RSS.*(* the model)
%             eval(['R_BIC_vesti(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').vestiPSTH3Dmodel.BIC_', models{m_inx}, ';']);
%             
%             % pack model fitting parameters to PARA.*.#
%             eval(['R_PARA_vesti(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').vestiPSTH3Dmodel.modelFitPara_', models{m_inx}, ';']);
%             
%         end
%     else
%         R_model(cell_inx).vestiSig = 0;
%         % Partial R_squared
%         R_PartR2_vesti(cell_inx).R2V = nan;
%         R_PartR2_vesti(cell_inx).R2A = nan;
%         R_PartR2_vesti(cell_inx).R2J = nan;
%         R_wVAJ_vesti(cell_inx).wV = nan;
%         R_wVAJ_vesti(cell_inx).wA = nan;
%         R_wVAJ_vesti(cell_inx).wJ = nan;
%         
%         for m_inx = 1:length(models)
%             % pack RSS values to RSS.*(* the model)
%             eval(['R_RSS_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack R_squared values to RSS.*(* the model)
%             eval(['R_Rsquared_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack BIC values to RSS.*(* the model)
%             eval(['R_BIC_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack model fitting parameters to PARA.*.#
%             eval(['R_PARA_vesti(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%         end
%     end
%     
%     if isstruct(R_model(cell_inx).visPSTH3Dmodel)
%         R_model(cell_inx).visSig = 1;
%         R_visNo = R_visNo+1;
%         % Partial R_squared
%         R_PartR2_vis(cell_inx).R2V = R_model(cell_inx).visPSTH3Dmodel.R2V;
%         R_PartR2_vis(cell_inx).R2A = R_model(cell_inx).visPSTH3Dmodel.R2A;
%         R_PartR2_vis(cell_inx).R2J = R_model(cell_inx).visPSTH3Dmodel.R2J;
%         R_wVAJ_vis(cell_inx).wV =  R_model(cell_inx).visPSTH3Dmodel.modelFitPara_VAJ(29)*(1-R_model(cell_inx).visPSTH3Dmodel.modelFitPara_VAJ(30));
%         R_wVAJ_vis(cell_inx).wA =  (1-R_model(cell_inx).visPSTH3Dmodel.modelFitPara_VAJ(29))*(1-R_model(cell_inx).visPSTH3Dmodel.modelFitPara_VAJ(30));
%         R_wVAJ_vis(cell_inx).wJ =  R_model(cell_inx).visPSTH3Dmodel.modelFitPara_VAJ(30);
%         
%         
%         for m_inx = 1:length(models)
%             % pack RSS values to RSS.*(* the model)
%             eval(['R_RSS_vis(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').visPSTH3Dmodel.rss_', models{m_inx} , ';']);
%             
%             % pack R_squared values to RSS.*(* the model)
%             eval(['R_Rsquared_vis(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').visPSTH3Dmodel.RSquared_', models{m_inx} , ';']);
%             
%             % pack BIC values to RSS.*(* the model)
%             eval(['R_BIC_vis(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').visPSTH3Dmodel.BIC_', models{m_inx}, ';']);
%             
%             % pack model fitting parameters to PARA.*.#
%             eval(['R_PARA_vis(',num2str(cell_inx),').', models{m_inx},' = R_model(',num2str(cell_inx),').visPSTH3Dmodel.modelFitPara_', models{m_inx}, ';']);
%             
%         end
%     else
%         R_model(cell_inx).visSig = 0;
%         % Partial R_squared
%         R_PartR2_vis(cell_inx).R2V = nan;
%         R_PartR2_vis(cell_inx).R2A = nan;
%         R_PartR2_vis(cell_inx).R2J = nan;
%         R_wVAJ_vis(cell_inx).wV = nan;
%         R_wVAJ_vis(cell_inx).wA = nan;
%         R_wVAJ_vis(cell_inx).wJ = nan;
%         
%         for m_inx = 1:length(models)
%             % pack RSS values to RSS.*(* the model)
%             eval(['R_RSS_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack R_squared values to RSS.*(* the model)
%             eval(['R_Rsquared_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack BIC values to RSS.*(* the model)
%             eval(['R_BIC_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%             % pack model fitting parameters to PARA.*.#
%             eval(['R_PARA_vis(',num2str(cell_inx),').', models{m_inx},' = nan;']);
%             
%         end
%     end
% end
% 
% disp('R model data loaded SUCCESS!');

% save the data
save('PSTH3DModel_T_OriData.mat','T_model','T_PARA_vesti', 'T_RSS_vesti', 'T_BIC_vesti','T_Rsquared_vesti','T_vestiNo','T_PartR2_vesti','T_wVAJ_vesti');
% save('PSTH3DModel_R_OriData.mat','R_model','R_PARA_vis', 'R_RSS_vis', 'R_BIC_vis','R_Rsquared_vis','R_PARA_vesti', 'R_RSS_vesti', 'R_BIC_vesti','R_Rsquared_vesti','R_vestiNo','R_visNo','R_PartR2_vis','R_PartR2_vesti','R_wVAJ_vis','R_wVAJ_vesti');
disp('DATA SAVED!');


%             w_paras = {{'wV'},...
%                 {'wA'},...
%                 {'wV','wA'},...
%                 {'wV','wJ'},...
%                 {'wA','wJ'},...
%                 {'wV','wA','wJ'},...
%                 {'wV','wA','wJ','wP'},...
%                 };
