% Population analysis for 3D model
% 101# Best model fit distribution according to BIC
% 102# R_squared distribution of each model
% 103# Partial R_squared distribution
% LBY 20171205

%% load data & pack data
% clear all;
cd('Z:\Data\TEMPO\BATCH\QQ_3DTuning');
load('Z:\Data\TEMPO\BATCH\QQ_3DTuning\PSTH3DModel_T_OriData.mat');
load('Z:\Data\TEMPO\BATCH\QQ_3DTuning\PSTH3DModel_R_OriData.mat');
Monkey = 'QQ';
models = {'VO','AO','VA','VJ','AJ','VAJ'};
%% analysis

colorDefsLBY;
T_vestiSig = cat(1,T_model.vestiSig);
T_visSig = cat(1,T_model.visSig);
R_vestiSig = cat(1,R_model.vestiSig);
R_visSig = cat(1,R_model.visSig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% basic infos.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp([Monkey,':']);
disp('Total cells for Translation 3D model: ');
disp(['[ vestibular: ',num2str(T_vestiNo),' ]   [  visual: ',num2str(T_visNo),' ]']);
disp('Total cells for Rotation 3D model: ');
disp(['[ vestibular: ',num2str(R_vestiNo),' ]   [  visual: ',num2str(R_visNo),' ]']);

%%%%%%%%%%%%%%%%%%%% Best model fit according to BIC %%%%%%%%%%%%%%%%%%%%%%

[~,T_BIC_min_inx_vesti] = min(squeeze(cell2mat(struct2cell(T_BIC_vesti))));
[~,T_BIC_min_inx_vis] = min(squeeze(cell2mat(struct2cell(T_BIC_vis))));
[~,R_BIC_min_inx_vesti] = min(squeeze(cell2mat(struct2cell(R_BIC_vesti))));
[~,R_BIC_min_inx_vis] = min(squeeze(cell2mat(struct2cell(R_BIC_vis))));

T_BIC_min_inx_vesti(T_vestiSig==0) = nan;
T_BIC_min_inx_vis(T_visSig==0) = nan;
R_BIC_min_inx_vesti(R_vestiSig==0) = nan;
R_BIC_min_inx_vis(R_visSig==0) = nan;

T_BIC_min_vesti_hist = hist(T_BIC_min_inx_vesti,length(models));
T_BIC_min_vis_hist = hist(T_BIC_min_inx_vis,length(models));
R_BIC_min_vesti_hist = hist(R_BIC_min_inx_vesti,length(models));
R_BIC_min_vis_hist = hist(R_BIC_min_inx_vis,length(models));

% % figures
% figure(101);set(gcf,'pos',[300 200 1000 600]);clf;
% BestFitModel = [T_BIC_min_vesti_hist;T_BIC_min_vis_hist;R_BIC_min_vesti_hist;R_BIC_min_vis_hist]';
% h = bar(BestFitModel,'grouped');
% xlabel('Models');ylabel('Cell Numbers (n)');
% set(gca,'xticklabel',models);
% set(h(1),'facecolor',colorDBlue,'edgecolor',colorDBlue);
% set(h(2),'facecolor',colorDRed,'edgecolor',colorDRed);
% set(h(3),'facecolor',colorLBlue,'edgecolor',colorLBlue);
% set(h(4),'facecolor',colorLRed,'edgecolor',colorLRed);
% title('Best model fit');
% h_l = legend(['Translation (Vestibular), n = ',num2str(T_vestiNo)],['Translation (Visual), n = ',num2str(T_visNo)],['Rotation (Vestibular), n = ',num2str(R_vestiNo)],['Rotation (Visual), n = ',num2str(R_visNo)],'location','NorthWest');
% set(h_l,'fontsize',15);
% SetFigure(15);
% set(gcf,'paperpositionmode','auto');
% saveas(101,'Z:\LBY\Population Results\BestFitModel','emf');

%%%%%%%%%%%%%%%%%  R_squared distribution of each model %%%%%%%%%%%%%%%%%%%

RSquared_T_vesti = squeeze(cell2mat(struct2cell(T_Rsquared_vesti)))';
RSquared_T_vis = squeeze(cell2mat(struct2cell(T_Rsquared_vis)))';
RSquared_R_vesti = squeeze(cell2mat(struct2cell(R_Rsquared_vesti)))';
RSquared_R_vis = squeeze(cell2mat(struct2cell(R_Rsquared_vis)))';

RSquared_T_vesti(T_vestiSig==0) = nan;
RSquared_T_vis(T_vestiSig==0) = nan;
RSquared_R_vesti(R_vestiSig==0) = nan;
RSquared_R_vis(R_vestiSig==0) = nan;

% figures
figure(102);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,3,0.1,0.15);

axes(h_subplot(1));
text(0.9,-0.3,'R^2','Fontsize',30,'rotation',90);
text(1.1,0.1,'Translation','Fontsize',25,'rotation',90);
text(1.1,-1.1,'Rotation','Fontsize',25,'rotation',90);
axis off;

axes(h_subplot(2));
hold on;
plot(RSquared_T_vesti','-o','color',colorLGray,'markeredgecolor',colorDBlue);
axis on;
xlim([0.5 6.5]);ylim([-0.5 1]);
set(gca,'xTick',1:6,'xticklabel',models);
title('Vestibular');

axes(h_subplot(3));
hold on;
plot(RSquared_T_vis','-o','color',colorLGray,'markeredgecolor',colorDRed);
axis on;
xlim([0.5 6.5]);ylim([-0.5 1]);
set(gca,'xTick',1:6,'xticklabel',models);
title('Visual');

axes(h_subplot(5));
hold on;
plot(RSquared_R_vesti','-o','color',colorLGray,'markeredgecolor',colorLBlue);
axis on;
xlim([0.5 6.5]);ylim([-0.5 1]);
set(gca,'xTick',1:6,'xticklabel',models);
xlabel('Models');

axes(h_subplot(6));
hold on;
plot(RSquared_R_vis','-o','color',colorLGray,'markeredgecolor',colorLRed);
axis on;
xlim([0.5 6.5]);ylim([-0.5 1]);
set(gca,'xTick',1:6,'xticklabel',models);
xlabel('Models');


SetFigure(20);
set(gcf,'paperpositionmode','auto');
saveas(102,'Z:\LBY\Population Results\RSquared_Distribution','emf');

%%%%%%%%%%%%%%%%%%%% Partial R_squared distribution  %%%%%%%%%%%%%%%%%%%%%%


PartR2_T_vesti = squeeze(cell2mat(struct2cell(T_PartR2_VAT_vesti)))';
PartR2_T_vis = squeeze(cell2mat(struct2cell(T_PartR2_VAT_vis)))';
PartR2_R_vesti = squeeze(cell2mat(struct2cell(R_PartR2_VAT_vesti)))';
PartR2_R_vis = squeeze(cell2mat(struct2cell(R_PartR2_VAT_vis)))';

PartR2_T_vesti(T_vestiSig==0) = nan;
PartR2_T_vis(T_vestiSig==0) = nan;
PartR2_R_vesti(R_vestiSig==0) = nan;
PartR2_R_vis(R_vestiSig==0) = nan;

% figures
figure(103);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,3,0.1,0.15);

axes(h_subplot(1));
text(0.9,-0.3,'Partial R^2','Fontsize',30,'rotation',90);
text(1.1,0.1,'Translation','Fontsize',25,'rotation',90);
text(1.1,-1.1,'Rotation','Fontsize',25,'rotation',90);
axis off;

axes(h_subplot(2));
hold on;
plot(PartR2_T_vesti','-o','color',colorLGray,'markeredgecolor',colorDBlue);
axis on;
xlim([0.5 3.5]);ylim([-0.5 1]);
set(gca,'xTick',1:3,'xticklabel',{'V/AJ','A/VJ','J/VA'});
title('Vestibular');

axes(h_subplot(3));
hold on;
plot(PartR2_T_vis','-o','color',colorLGray,'markeredgecolor',colorDRed);
axis on;
xlim([0.5 3.5]);ylim([-0.5 1]);
set(gca,'xTick',1:3,'xticklabel',{'V/AJ','A/VJ','J/VA'});
title('Visual');

axes(h_subplot(5));
hold on;
plot(PartR2_R_vesti','-o','color',colorLGray,'markeredgecolor',colorLBlue);
axis on;
xlim([0.5 3.5]);ylim([-0.5 1]);
set(gca,'xTick',1:3,'xticklabel',{'V/AJ','A/VJ','J/VA'});
xlabel('Models');

axes(h_subplot(6));
hold on;
plot(PartR2_R_vis','-o','color',colorLGray,'markeredgecolor',colorLRed);
axis on;
xlim([0.5 3.5]);ylim([-0.5 1]);
set(gca,'xTick',1:3,'xticklabel',{'V/AJ','A/VJ','J/VA'});
xlabel('Models');

SetFigure(20);
set(gcf,'paperpositionmode','auto');
saveas(103,'Z:\LBY\Population Results\Partial_RSquared_Distribution','emf');

%%%%%%%%%%%%%%%%%%%  R_squared distribution  %%%%%%%%%%%%%%%%%%%%%%

models = {'VO','AO','VA','VJ','AJ','VAJ'};
xR2 = 0.05:0.1:0.75;

% figures
figure(104);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,4,0.1,0.15);

for ii = 1:6
    
axes(h_subplot(ii));
hold on;
[nelements, ncenters] = hist(RSquared_T_vesti(:,ii),xR2);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(T_vestiSPeakT_plot))]);
plot(nanmedian(RSquared_T_vesti(:,ii)),max(nelements)*1.1,'kv');
text(nanmedian(RSquared_T_vesti(:,ii))*1.1,max(nelements)*1.2,num2str(nanmedian(RSquared_T_vesti(:,ii))));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Single-peaked');
title([models{ii},' model']);
axis on;
hold off;
          
end

suptitle('Translation - vestibular');
SetFigure(15);

figure(105);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,4,0.1,0.15);

for ii = 1:6
    
axes(h_subplot(ii));
hold on;
[nelements, ncenters] = hist(RSquared_T_vis(:,ii),xR2);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(T_visSPeakT_plot))]);
plot(nanmedian(RSquared_T_vis(:,ii)),max(nelements)*1.1,'kv');
text(nanmedian(RSquared_T_vis(:,ii))*1.1,max(nelements)*1.2,num2str(nanmedian(RSquared_T_vis(:,ii))));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Single-peaked');
title([models{ii},' model']);
axis on;
hold off;
          
end

suptitle('Translation - visual');
SetFigure(15);


% figures
figure(106);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,4,0.1,0.15);

for ii = 1:6
    
axes(h_subplot(ii));
hold on;
[nelements, ncenters] = hist(RSquared_R_vesti(:,ii),xR2);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(R_vestiSPeakR_plot))]);
plot(nanmedian(RSquared_R_vesti(:,ii)),max(nelements)*1.1,'kv');
text(nanmedian(RSquared_R_vesti(:,ii))*1.1,max(nelements)*1.2,num2str(nanmedian(RSquared_R_vesti(:,ii))));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Single-peaked');
title([models{ii},' model']);
axis on;
hold off;
          
end

suptitle('Rotation - vestibular');
SetFigure(15);

figure(107);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,4,0.1,0.15);

for ii = 1:6
    
axes(h_subplot(ii));
hold on;
[nelements, ncenters] = hist(RSquared_R_vis(:,ii),xR2);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(R_visSPeakR_plot))]);
plot(nanmedian(RSquared_R_vis(:,ii)),max(nelements)*1.1,'kv');
text(nanmedian(RSquared_R_vis(:,ii))*1.1,max(nelements)*1.2,num2str(nanmedian(RSquared_R_vis(:,ii))));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Single-peaked');
title([models{ii},' model']);
axis on;
hold off;
          
end

suptitle('Rotation - visual');
SetFigure(15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weight for VAJ model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(105);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,3,[0.1 0.02],0.15);

axes(h_subplot(2));
T_vesti_w = squeeze(cell2mat(struct2cell(T_wVAJ_vesti)))';
%-- Plot the axis system
[h,hg,htick]=terplot(5);
%-- Plot the data ...
hter=ternaryc(T_vesti_w(:,3),T_vesti_w(:,1),T_vesti_w(:,2));
%-- ... and modify the symbol:
set(hter,'marker','o','markerfacecolor',colorDBlue,'markersize',7,'markeredgecolor','w');
terlabel('wJ','wV','wA');
% view(180,-90);
title('T - vestibular');
% axis off;

axes(h_subplot(3));
T_vis_w = squeeze(cell2mat(struct2cell(T_wVAJ_vis)))';
%-- Plot the axis system
[h,hg,htick]=terplot(5);
%-- Plot the data ...
hter=ternaryc(T_vis_w(:,3),T_vis_w(:,1),T_vis_w(:,2));
%-- ... and modify the symbol:
set(hter,'marker','o','markerfacecolor',colorDRed,'markersize',7,'markeredgecolor','w');
terlabel('wJ','wV','wA');

title('T - visual');
% axis off;

axes(h_subplot(5));
R_vesti_w = squeeze(cell2mat(struct2cell(R_wVAJ_vesti)))';
%-- Plot the axis system
[h,hg,htick]=terplot(5);
%-- Plot the data ...
hter=ternaryc(R_vesti_w(:,3),R_vesti_w(:,1),R_vesti_w(:,2));
%-- ... and modify the symbol:
set(hter,'marker','o','markerfacecolor',colorLBlue,'markersize',7,'markeredgecolor','w');
terlabel('wJ','wV','wA');

title('R - vestibular');
% axis off;

axes(h_subplot(6));
R_vis_w = squeeze(cell2mat(struct2cell(R_wVAJ_vis)))';
%-- Plot the axis system
[h,hg,htick]=terplot(5);
%-- Plot the data ...
hter=ternaryc(R_vis_w(:,3),R_vis_w(:,1),R_vis_w(:,2));
%-- ... and modify the symbol:
set(hter,'marker','o','markerfacecolor',colorLRed,'markersize',7,'markeredgecolor','w');
terlabel('wJ','wV','wA');

title('R - visual');
% axis off;

SetFigure(15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% weight for VA model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(106);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,3,[0.1 0.02],0.15);

axes(h_subplot(2));hold on;
T_vesti_w = squeeze(cell2mat(struct2cell(T_wVA_vesti)))';
plot(T_vesti_w(:,1),T_vesti_w(:,2),'o','markerfacecolor',colorDBlue,'markersize',7,'markeredgecolor','w');
xlabel('wV');ylabel('wA');
plot([0 1],[0.5 0.5],'--','color',colorLGray);
plot([0.5 0.5],[0 1],'--','color',colorLGray);
axis square;
% view(180,-90);
title('T - vestibular');
axis on;

axes(h_subplot(3));hold on;
T_vis_w = squeeze(cell2mat(struct2cell(T_wVA_vis)))';
plot(T_vis_w(:,1),T_vis_w(:,2),'o','markerfacecolor',colorDRed,'markersize',7,'markeredgecolor','w');
plot([0 1],[0.5 0.5],'--','color',colorLGray);
plot([0.5 0.5],[0 1],'--','color',colorLGray);
xlabel('wV');ylabel('wA');
axis square;
% view(180,-90);
title('T - visual');
axis on;

axes(h_subplot(5));hold on;
R_vesti_w = squeeze(cell2mat(struct2cell(R_wVA_vesti)))';
plot(R_vesti_w(:,1),R_vesti_w(:,2),'o','markerfacecolor',colorLBlue,'markersize',7,'markeredgecolor','w');
xlabel('wV');ylabel('wA');
plot([0 1],[0.5 0.5],'--','color',colorLGray);
plot([0.5 0.5],[0 1],'--','color',colorLGray);
axis square;
% view(180,-90);
title('R - vestibular');
axis on;

axes(h_subplot(6));hold on;
R_vis_w = squeeze(cell2mat(struct2cell(R_wVA_vis)))';
plot(R_vis_w(:,1),R_vis_w(:,2),'o','markerfacecolor',colorLRed,'markersize',7,'markeredgecolor','w');
plot([0 1],[0.5 0.5],'--','color',colorLGray);
plot([0.5 0.5],[0 1],'--','color',colorLGray);
xlabel('wV');ylabel('wA');
axis square;
% view(180,-90);
title('R - visual');
axis on;

SetFigure(15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% weight for VA model (ratio distribution) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xRatio = -2:0.4:2;

figure(107);set(gcf,'pos',[60 70 1500 800]);clf;
[~,h_subplot] = tight_subplot(2,3,[0.2 0.1],0.15);

axes(h_subplot(2));hold on;
T_vesti_w = squeeze(cell2mat(struct2cell(T_wVA_vesti)))';
T_vesti_ratio = log(T_vesti_w(:,1)./T_vesti_w(:,2));
[nelements, ncenters] = hist(T_vesti_ratio,xRatio);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(T_visDPeakT_plot))]);
plot(nanmedian(T_vesti_ratio),max(nelements)*1.1,'kv');
text(nanmedian(T_vesti_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(T_vesti_ratio)));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Double-peaked, early');
% [p,h] = ranksum(T_vesti_ratio,1)
title('T-vestibular');
xlim([-2 2]);
% set(gca,'xscal','log');
axis on;
hold off;

axes(h_subplot(3));hold on;
T_vis_w = squeeze(cell2mat(struct2cell(T_wVA_vis)))';
T_vis_ratio = log(T_vis_w(:,1)./T_vis_w(:,2));
[nelements, ncenters] = hist(T_vis_ratio,xRatio);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(T_visDPeakT_plot))]);
plot(nanmedian(T_vis_ratio),max(nelements)*1.1,'kv');
text(nanmedian(T_vis_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(T_vis_ratio)));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Double-peaked, early');
axis on;
hold off;
title('T-visual');
xlim([-2 2]);
[p,h] = ranksum(T_vis_ratio,1)

axes(h_subplot(5));hold on;
R_vesti_w = squeeze(cell2mat(struct2cell(R_wVA_vesti)))';
R_vesti_ratio = log(R_vesti_w(:,1)./R_vesti_w(:,2));
[nelements, ncenters] = hist(R_vesti_ratio,xRatio);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(R_visDPeakR_plot))]);
plot(nanmedian(R_vesti_ratio),max(nelements)*1.1,'kv');
text(nanmedian(R_vesti_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(R_vesti_ratio)));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Double-peaked, early');
axis on;
hold off;
title('R-vestibular');
xlim([-2 2]);
[p,h] = ranksum(R_vesti_ratio,1)

axes(h_subplot(6));hold on;
R_vis_w = squeeze(cell2mat(struct2cell(R_wVA_vis)))';
R_vis_ratio = log(R_vis_w(:,1)./R_vis_w(:,2));
[nelements, ncenters] = hist(R_vis_ratio,xRatio);
h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
set(h1,'linewidth',1.5);
% text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(R_visDPeakR_plot))]);
plot(nanmedian(R_vis_ratio),max(nelements)*1.1,'kv');
text(nanmedian(R_vis_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(R_vis_ratio)));
% set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% xlabel('Double-peaked, early');
axis on;
hold off;
title('R-visual');
xlim([-2 2]);
[p,h] = ranksum(R_vis_ratio,1)

SetFigure(15);

% %%%%%%%%%%%%% compared with other areas %%%%%%%%%%%%%%%%%%%%%%%%
% % data from Laurens, VAJ model, normalized
% 
% xRatio = -4:0.5:4;
% 
% figure(108);set(gcf,'pos',[60 70 1500 800]);clf;
% [~,h_subplot] = tight_subplot(2,5,[0.2 0.02],0.15);
% 
% axes(h_subplot(1));hold on;
% OA_ratio = log(OA_N(:,1)./OA_N(:,2));
% [nelements, ncenters] = hist(OA_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(OA_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(OA_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(OA_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% [p,h] = ranksum(OA_ratio,0)
% title('OA');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(3));hold on;
% VN_ratio = log(VN_N(:,1)./VN_N(:,2));
% [nelements, ncenters] = hist(VN_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(VN_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(VN_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(VN_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(VN_ratio,1)
% title('VN');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(4));hold on;
% CN_ratio = log(CN_N(:,1)./CN_N(:,2));
% [nelements, ncenters] = hist(CN_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(CN_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(CN_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(CN_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(CN_ratio,1)
% title('CN');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(5));hold on;
% T_vesti_w = squeeze(cell2mat(struct2cell(T_wVA_vesti)))';
% T_vesti_ratio = log(T_vesti_w(:,1)./T_vesti_w(:,2));
% [nelements, ncenters] = hist(T_vesti_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% % text(170,max(max(nelements),max(nelements)),['n = ',num2str(length(T_visDPeakT_plot))]);
% plot(nanmedian(T_vesti_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(T_vesti_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(T_vesti_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(T_vesti_ratio,1)
% title('PCC');
% % % xlim([-1 5]);
% % set(gca,'xscal','log');
% axis on;
% hold off;
% 
% 
% axes(h_subplot(6));hold on;
% PIVC_ratio = log(PIVC_N(:,1)./PIVC_N(:,2));
% [nelements, ncenters] = hist(PIVC_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(PIVC_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(PIVC_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(PIVC_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(PIVC_ratio,1)
% title('PIVC');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(7));hold on;
% VPS_ratio = log(VPS_N(:,1)./VPS_N(:,2));
% [nelements, ncenters] = hist(VPS_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(VPS_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(VPS_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(VPS_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(VPS_ratio,1)
% title('VPS');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(8));hold on;
% MST_ratio = log(MST_N(:,1)./MST_N(:,2));
% [nelements, ncenters] = hist(MST_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(MST_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(MST_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(MST_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% [p,h] = ranksum(MST_ratio,1)
% title('MST');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(9));hold on;
% VIP_ratio = log(VIP_N(:,1)./VIP_N(:,2));
% [nelements, ncenters] = hist(VIP_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(VIP_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(VIP_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(VIP_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(VIP_ratio,1)
% title('VIP');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% axes(h_subplot(10));hold on;
% FEF_ratio = log(FEF_N(:,1)./FEF_N(:,2));
% [nelements, ncenters] = hist(FEF_ratio,xRatio);
% h1 = bar(ncenters, nelements, 0.7,'k','edgecolor','k');
% set(h1,'linewidth',1.5);
% plot(nanmedian(FEF_ratio),max(nelements)*1.1,'kv');
% text(nanmedian(FEF_ratio)*1.1,max(nelements)*1.2,num2str(nanmedian(FEF_ratio)));
% % set(gca,'xtick',[0 500 1000 1500],'xticklabel',[],'xlim',[0 1600]);
% % xlabel('Double-peaked, early');
% % [p,h] = ranksum(FEF_ratio,1)
% title('FEF');
% % xlim([-1 5]);
% axis on;
% hold off;
% 
% SetFigure(15);


% %%%%%%%%%%%%%%%%%%%%%%%%% Weight vs. R2 (VA model) %%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(108);set(gcf,'pos',[60 70 1200 800]);clf;
% [~,h_subplot] = tight_subplot(2,3,[0.2 0.2],0.15);
% 
% axes(h_subplot(2));hold on;
% T_vesti_w = squeeze(cell2mat(struct2cell(T_wVA_vesti)))';
% T_vesti_ratio = log(T_vesti_w(:,1)./T_vesti_w(:,2));
% T_vesti_ratio = T_vesti_ratio(~isnan(T_vesti_ratio));
% T_vesti_r2 = cat(1,T_Rsquared_vesti.VA);
% T_vesti_r2 = T_vesti_r2(~isnan(T_vesti_r2));
% plot(T_vesti_ratio,T_vesti_r2,'o','markersize',6,'markerfacecolor',colorDBlue);
% T_vesti = polyfit(T_vesti_ratio, T_vesti_r2, 1);
% plot(T_vesti_ratio, polyval(T_vesti, T_vesti_ratio),'color',colorDBlue,'linewidth',3);
% [T_r_vesti,T_p_vesti] =corrcoef(T_vesti_ratio,T_vesti_r2);
% title(sprintf('T-vestibular, r = %g \np = %g',T_r_vesti(1,2),T_p_vesti(1,2)));
% axis on;
% ylabel('R^2');
% % title('T-vestibular');
% 
% axes(h_subplot(3));hold on;
% T_vis_w = squeeze(cell2mat(struct2cell(T_wVA_vis)))';
% T_vis_ratio = log(T_vis_w(:,1)./T_vis_w(:,2));
% T_vis_ratio = T_vis_ratio(~isnan(T_vis_ratio));
% T_vis_r2 = cat(1,T_Rsquared_vis.VA);
% T_vis_r2 = T_vis_r2(~isnan(T_vis_r2));
% plot(T_vis_ratio,T_vis_r2,'o','markersize',6,'markerfacecolor',colorDRed);
% T_vis = polyfit(T_vis_ratio, T_vis_r2, 1);
% plot(T_vis_ratio, polyval(T_vis, T_vis_ratio),'color',colorDRed,'linewidth',3);
% [T_r_vis,T_p_vis] =corrcoef(T_vis_ratio,T_vis_r2);
% title(sprintf('T-visual, r = %g \np = %g',T_r_vis(1,2),T_p_vis(1,2)));
% axis on;
% ylabel('R^2');
% % title('T-visual');
% 
% 
% 
% axes(h_subplot(5));hold on;
% R_vesti_w = squeeze(cell2mat(struct2cell(R_wVA_vesti)))';
% R_vesti_ratio = log(R_vesti_w(:,1)./R_vesti_w(:,2));
% R_vesti_ratio = R_vesti_ratio(~isnan(R_vesti_ratio));
% R_vesti_r2 = cat(1,R_Rsquared_vesti.VA);
% R_vesti_r2 = R_vesti_r2(~isnan(R_vesti_r2));
% plot(R_vesti_ratio,R_vesti_r2,'o','markersize',6,'markerfacecolor',colorLBlue);
% R_vesti = polyfit(R_vesti_ratio, R_vesti_r2, 1);
% plot(R_vesti_ratio, polyval(R_vesti, R_vesti_ratio),'color',colorLBlue,'linewidth',3);
% [R_r_vesti,R_p_vesti] =corrcoef(R_vesti_ratio,R_vesti_r2);
% title(sprintf('R-vestibular,r = %g \np = %g',R_r_vesti(1,2),R_p_vesti(1,2)));
% axis on;
% ylabel('R^2');
% 
% 
% axes(h_subplot(6));hold on;
% R_vis_w = squeeze(cell2mat(struct2cell(R_wVA_vis)))';
% R_vis_ratio = log(R_vis_w(:,1)./R_vis_w(:,2));
% R_vis_ratio = R_vis_ratio(~isnan(R_vis_ratio));
% R_vis_r2 = cat(1,R_Rsquared_vis.VA);
% R_vis_r2 = R_vis_r2(~isnan(R_vis_r2));
% plot(R_vis_ratio,R_vis_r2,'o','markersize',6,'markerfacecolor',colorLRed);
% R_vis = polyfit(R_vis_ratio, R_vis_r2, 1);
% plot(R_vis_ratio, polyval(R_vis, R_vis_ratio),'color',colorLRed,'linewidth',3);
% [R_r_vis,R_p_vis] =corrcoef(R_vis_ratio,R_vis_r2);
% title(sprintf('R-visual,r = %g \np = %g',R_r_vis(1,2),R_p_vis(1,2)));
% axis on;
% ylabel('R^2');
% % title('T-visual');
% 
% 
% SetFigure(15);
% 
% 
