%%   Microstimulation Pooled Data
% HH 201406
function Group_MicroStim(XlsData)

% %% Read xls
% [num,txt,raw] = xlsread('Result.xlsm',2);
% 
% % Get Header infomation
% HEADS_N = 3;
% 
% header_all = txt(HEADS_N-1,:);
% header_all(strcmp(header_all,'')) = txt(HEADS_N-2,strcmp(header_all,''));
% 
% for i = 1:length(header_all)
%     try
%         if i == num(1,i)
%             eval([header_all{i} '=' num2str(i) ';']);
%         else
%             disp('Header info error...');
%             keyboard;
%         end
%     catch
%     end
% end
% 
% % Delete headers
% end_line = find(~isnan(num(:,1)),1,'last'); 
% num = num(HEADS_N+1:end_line,:);
% txt = txt(HEADS_N : end_line - 1,:); % Note here
% raw = raw(HEADS_N+1:end_line,:);


% psth = XlsCell2Mat(raw(HEADS_N+1:end,PSTH));
% psth_t = str2num(raw{3,4});
% cp_shift = XlsCell2Mat(raw(HEADS_N+1:end,CP));
% nt_shift = XlsCell2Mat(raw(HEADS_N+1:end,NT));
% cp_t = str2num(raw{end,20});
% nt_t = cp_t;

num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;

%% Get data

% Mask
mask_all = strcmp(txt(:,header.Protocol),'u-stim');

ustim_num = num(mask_all,:);
ustim_txt = txt(mask_all,:); % Note here
ustim_raw = raw(mask_all,:);


% Areas and stim_types of interest
area_list = {'VIP','MST'};
area_list_color = {'b','r'};
stim_type_list = {'Vestibular','Vis'};  % Don't change the order: Vesti, Vis, Comb
stim_type_list_color = {'b','r','g'};

% -------------  Get data

% Basics
areas = ustim_txt(:,header.Area);
amps = ustim_num(:,header.uA);
reps = ustim_num(:,header.rep_ustim);

dPSE = ustim_num(:,[header.dPSE_vest header.dPSE_vis header.dPSE_comb]);
rSIG = ustim_num(:,[header.rSigma_vest  header.rSigma_vis  header.rSigma_comb]);

pPSE = ustim_num(:,[header.PSE_p_vest header.PSE_p_vis header.PSE_p_comb]);
pSIG = ustim_num(:,[header.sigma_p_vest header.sigma_p_vis header.sigma_p_comb]);

tuningOK = ustim_num(:,[header.tuning_ok_vest header.tuning_ok_vis header.tuning_ok_comb]);
electrodeOK = ustim_num(:,header.electrode_failure);

% Tuning indices
tuningDDI = ustim_num(:,[header.DDI_vest header.DDI_vis header.DDI_comb]);
tuningHTI = ustim_num(:,[header.HTI_vest header.HTI_vis header.HTI_comb]);
tuningdPrime = ustim_num(:,[header.dPrime0_vest header.dPrime0_vis header.dPrime0_comb]);

% Tuning alignment  HH20141216
tuningPref = ustim_num(:,[header.Pref_vest header.Pref_vis header.Pref_comb]);

% Mean tuning curves
tuning0 = XlsCell2Mat(ustim_raw(:,[header.mean_vest header.mean_vis header.mean_comb]));
tuning1 = XlsCell2Mat(ustim_raw(:,[header.mean_vest1 header.mean_vis1 header.mean_comb1]));
tuning2 = XlsCell2Mat(ustim_raw(:,[header.mean_vest2 header.mean_vis2 header.mean_comb2]));

tuningClustering = zeros(size(tuningDDI))*NaN;

fig_number = 1;

for a = 1:length(area_list)
    mask_a = strcmp(areas,area_list{a});
    
    for k = 1:length(stim_type_list)
        
        mask_a_k = mask_a & ~isnan(dPSE(:,k));
        

        %% PSE shift

        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k} ', dPSE'], 'Position',[300 300 410 300]); clf
        fig_number = fig_number + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % delta_PSE MASK
        
        mask_PSE = mask_a_k & (0|reps >= 8) & (0|tuningOK(:,k) > 0) & (0|amps <= 20) & (0 | electrodeOK == 0 | electrodeOK ==3);
        
        if sum(mask_PSE) == 0 ; continue; end
        
        S_PSE = mask_PSE & (pPSE(:,k) < 0.05);
        NS_PSE= mask_PSE & (pPSE(:,k) >= 0.05);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Flip dPSE if pref is not aligned with Combine direction  HH20141216
        toFlip =  (sign(tuningPref(mask_PSE,1)) ~= sign(tuningPref(mask_PSE,3))) & ~isnan(tuningPref(mask_PSE,3));
        find_PSE = find(mask_PSE);
%         dPSE(find_PSE(toFlip)) = -dPSE(find_PSE(toFlip));
                
        bin = 0.5;
        xcenters = floor(min(dPSE(S_PSE|NS_PSE,k))/2)*2-bin/2:bin:ceil(max(dPSE(S_PSE|NS_PSE,k))/2)*2;  % Ensure that 0 is one of the borders.
        
        histS = hist(dPSE(S_PSE,k),xcenters);
        histNS = hist(dPSE(NS_PSE,k),xcenters);
        
        hbars = bar(xcenters,[histS' histNS'],1,'stacked','k','LineWidth',2);
        set(hbars,'EdgeColor',stim_type_list_color{k},'FaceColor',stim_type_list_color{k});
        set(hbars(2),'FaceColor','none');
        
        
        all_PSE = dPSE(S_PSE|NS_PSE,k);
        meanPSE = mean(all_PSE);
        maxY = max(max([histS+histNS])*1.2,6);
        axis([min(xcenters)-1 max(xcenters)+1 0 maxY]);
        set(gca,'xtick',[-100:100],'xMinorTick','on');
        xlim([min(-3,min(all_PSE))-0.5 max(3,max(all_PSE))+0.5]);
        
        % ttest
        [~,p] = ttest(all_PSE,0);
        
        % Annotation        
        hold on; 
        plot([meanPSE meanPSE],[maxY/1 maxY/1.2],stim_type_list_color{k}); 
        fill([meanPSE meanPSE+0.2 meanPSE-0.2],[maxY/1.2 maxY/1.2+0.5 maxY/1.2+0.5],stim_type_list_color{k},'EdgeColor',stim_type_list_color{k});
        tt = text(meanPSE+1, maxY/1.2+0.5, sprintf('Mean \\Delta\\mu = %2.2f^\\circ\n\\itp \\rm= %.3g\n\\Delta\\mu > 0: %g / %g (%.3g%%)\n\\itp \\rm< 0.05: %g / %g (%.3g%%)',...
            meanPSE,p, sum(all_PSE >0), length(all_PSE), sum(all_PSE>0)/length(all_PSE)*100, sum(S_PSE),sum(S_PSE)+sum(NS_PSE),sum(S_PSE)/(sum(S_PSE)+sum(NS_PSE))*100));
        set(tt,'Color',stim_type_list_color{k});
        
        plot([0 0],[0 maxY],'k:' );
        
        xlabel(sprintf('\\Delta\\mu (\\circ)'));
        ylabel('Number of cases');
        legend ({'p<0.05' 'p>0.05'},'location','best')
        SetFigure(15)
        
        
% %{
        %% 2-D Plot with different conditions
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2_D Plot MASK
        
        % all conditions
        mask_2D = mask_a_k &   (1|amps <= 50) & (1|reps >= 8) & (1|tuningOK(:,k) > 0) & (electrodeOK == 0 | electrodeOK ==3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k} ', 2-D Plot']); clf
        fig_number = fig_number + 1;
        
        result_200 = [];
        result_50_100 = [];
        result_20_20 = [];
        result_20_45 = [];
        
        for i = find(mask_2D)'
            
            if amps(i) <= 50
                %         if ~(tuningOK(i)==1); continue; end;
                markerSize = 10;
%             elseif  amps(i) <= 100
%                 markerSize = 20;
%                 result_50_100 = [result_50_100; dPSE(i,k) rSIG(i,k)];
            elseif amps(i) > 50
                markerSize = 40;
                result_200 = [result_200; dPSE(i,k) rSIG(i,k)];
            end
            
            if pPSE(i,k) < 0.05 && pSIG(i,k) < 0.05
                markerCol = [1 0 1];
                %         markerFaceCol = 'm';
            elseif pPSE(i,k) <= 0.05 && pSIG(i,k) >= 0.05
                markerCol = [1 0 0];
                %         markerFaceCol = 'r';
            elseif pPSE(i,k) >= 0.05 && pSIG(i,k) < 0.05
                markerCol =  [0 0 1];
                %         markerFaceCol = 'b';
            else
                markerCol = [.6 .6 .6];%'k';
                %         markerFaceCol = 'none';
            end
            
            markerFaceCol = 'none';
            marker = 'o';
            
            plot(dPSE(i,k),rSIG(i,k),[marker],'color',markerCol,'markersize',markerSize,'linewidth',2,'markerFaceCol',markerFaceCol);
            hold on;
            
        end
        
        set(gca,'yScale','log','ytick',[0.01 0.5 1 2 3 4 5],'yticklabel','0.01| 0.5| 1|||| 5');
        
        axis([-6 16 0.6 5.5]);
        
        lims = axis;
%         axis([min(-4,lims(1)) max(10,lims(2)) min(0.3,lims(3)) max(5,lims(4)) ]);
        lims = axis;
        
        text(lims(2)*0.8, lims(4)*0.8, sprintf('\\itN\\rm = %g',sum(mask_2D)));
        
        plot([lims(1)*1.1 lims(2)*1.1],[1 1],'k:');
        plot([0 0],[lims(3)/1.2 lims(4)*1.2],'k:');
        axis([lims(1:2)*1.1 lims(3)/1.2 lims(4)*1.2]);
        SetFigure(20);
        xlabel(sprintf('\\Delta\\mu (\\circ)'),'color','r');
        ylabel(sprintf('\\sigma_{ratio}'),'color','b');
        
        % axis([-6 15 0.5 6]);

%}
        
% %{
        %% dPSE v.s. Tuning properties
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dPSE v.s. Tuning properties MASK
        
        % 20140512 MST_100uA
        mask_Discri = mask_a_k & (0|amps <= 20) & (0|reps >= 8) & (0|tuningOK(:,k) > 0) & (electrodeOK == 0 | electrodeOK ==3);
        
        % 20140512 VIP_20uA, d'0degree
        % mask_Discri = (0|tuningOK > 0.5) & (0|amp <= 50);
        
        AMP_20_S = mask_Discri & (pPSE(:,k)<0.05);
        AMP_20_NS = mask_Discri & (pPSE(:,k)>=0.05);
        % AMP_20_S = AMP_20_S(1:size(tuning,1));
        % AMP_20_NS = AMP_20_NS(1:size(tuning,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        % ----------------- 1. DDI/HTI and dPSE (linear)-----------------
        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k}]); clf
        fig_number = fig_number + 1;

        % DDI
        DDI_S = abs(tuningDDI(AMP_20_S,k));
        DDI_NS = abs(tuningDDI(AMP_20_NS,k));
        
        % Deviation of pref. heading from +/- 90
        % dis_S = abs((abs(tuning(AMP_20_S,1))-90));
        % dis_NS = abs((abs(tuning(AMP_20_NS,1))-90));
        
        
        plot(DDI_S,abs(dPSE(AMP_20_S,k)),'ro','markerfacecolor','r');
        hold on;  plot(DDI_NS,abs(dPSE(AMP_20_NS,k)),'ro');
        set(findall(gca,'color','r'),'markersize',10);
        set(findall(gca,'color','r'),'linewidth',2);
        set(gca,'ticklength',[0.02 0]);
        
        % fitting
        allX = ([DDI_S; DDI_NS]);
        allY = (abs([dPSE(AMP_20_S,k); dPSE(AMP_20_NS,k)]));
        
        axis([min(allX)/1.5 max(allX)*1.5 min(allY)/1.5 max(allY)*1.5]);
        
        [r,p]=corrcoef(allX,allY);
        [para,S]=polyfit(allX,allY,1);
        xx = min(allX)-0.5:0.1:max(allX)+0.5;
        Y = polyval(para,xx);
        plot(xx,Y,'r-','linewidth',2);
        
        xlabel('DDI');
        ylabel(sprintf('|\\Delta\\mu (\\circ)|'));
        title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))]);
        SetFigure(20);

        % ----------------- 2.D-prime and dPSE (log scale) -----------------
        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k} ', dPrime']); clf
        fig_number = fig_number + 1;
        %
        % dis_S = (abs(tuning(AMP_20_S,6))+abs(tuning(AMP_20_S,5)))/2;
        % dis_NS = (abs(tuning(AMP_20_NS,6))+abs(tuning(AMP_20_NS,5)))/2;
        %
        dPrime_S = abs(tuningdPrime(AMP_20_S,k));
        dPrime_NS = abs(tuningdPrime(AMP_20_NS,k));
        
%         loglog(dPrime_S,abs(dPSE(AMP_20_S,k)),'ko','markerfacecolor','k');

        plot(dPrime_S,dPSE(AMP_20_S,k),'ko','markerfacecolor','k'); set(gca,'xscale','log');
        hold on;  plot(dPrime_NS,dPSE(AMP_20_NS,k),'ko');
        set(findall(gca,'color','k'),'markersize',10);
        set(findall(gca,'color','k'),'linewidth',2);
        set(gca,'ticklength',[0.02 0]);
        
        % fitting
        allX = ([dPrime_S; dPrime_NS]);
        allY = ([dPSE(AMP_20_S,k); dPSE(AMP_20_NS,k)]);
        
        [r,p]=corrcoef(log(allX),allY);
        [para,S]=polyfit(log(allX),allY,1);
        xx = min(log(allX))-0.5:0.1:max(log(allX))+0.5;
        Y = polyval(para,xx);
        plot(exp(xx),Y,'k-','linewidth',2);
        
        axis([min(allX)/1.5 max(allX)*1.5 min(allY)/1.5 max(allY)*1.5]);
        xlabel('Discriminability (d'')');
        ylabel(sprintf('\\Delta\\mu (\\circ)'));
        title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))]);
        set(gca,'xtick',[0.1 :0.1 1:1:10]);
        SetFigure(20);
        
        LinearCorrelation({dPrime_S,dPrime_NS},{dPSE(AMP_20_S,k),dPSE(AMP_20_NS,k)},'logx',1,...
        'Xlabel','Discriminability (d'')','Ylabel',sprintf('\\Delta\\mu (\\circ)'),...
        'FaceColors',{'k','none'},'Markers',{'o'},'Method','Pearson','FittingMethod',2, ...
        'LineStyles',{'k-','k:'},'MarkerSize',12,'CombinedIndex',[3],'figN',fig_number); fig_number = fig_number + 1;   


        % ----------------- 3.Local clustering -----------------
        for i = find(AMP_20_S | AMP_20_NS)'
            tuning_this{1} = tuning0{k}(i,:);   % 0 um
            tuning_this{2} = tuning1{k}(i,:);   % -100 um
            tuning_this{3} = tuning2{k}(i,:);   % 100 um
            
            % Use 45 degree tuning curve to calculate the correlation
            for j = 1:3
                len_azi = sum(~isnan(tuning_this{j}));
                switch len_azi
                    case 8
                        ind = 1:8;
                    case 10
                        ind = setdiff(1:10,[3 5]);
                    case 16
                        ind = 1:2:16;
                    case 0
                        ind = [];
                    otherwise
                        disp('Wrong tuning length...');
                        keyboard
                end
                tuning_this{j} = tuning_this{j}(ind);
            end
            
            if ~isempty(tuning_this{2})
                [r1,~] = corrcoef(tuning_this{1},tuning_this{2});
            else
                r1 = [NaN NaN ; NaN NaN];
            end
            
            
            if ~isempty(tuning_this{3})
                [r2,~] = corrcoef(tuning_this{1},tuning_this{3});
            else
                r2 = r1;
            end
            
            tuningClustering(i,k) = (r1(2) + r2(2))/2;
        end
        
        % ----  Plot Clustering ---
        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k} ', Clustering']); clf
        fig_number = fig_number + 1;

        % Local clustering
        LocCluster_S = tuningClustering(AMP_20_S,k);
        LocCluster_NS = tuningClustering(AMP_20_NS,k);
                
        plot(LocCluster_S,abs(dPSE(AMP_20_S,k)),'ko','markerfacecolor','k');
        hold on;  plot(LocCluster_NS,abs(dPSE(AMP_20_NS,k)),'ko');
        set(findall(gca,'color','k'),'markersize',10);
        set(findall(gca,'color','k'),'linewidth',2);
        set(gca,'ticklength',[0.02 0]);
        
        % fitting
        allX = ([LocCluster_S; LocCluster_NS]);
        allY = (abs([dPSE(AMP_20_S,k); dPSE(AMP_20_NS,k)]));
        
        axis auto;
        % axis([0 1 min(allY)/1.5 max(allY)*1.5]);
        
        [r,p]=corrcoef(allX,allY);
        [para,S]=polyfit(allX,allY,1);
        xx = min(allX)-0.5:0.1:max(allX)+0.5;
        Y = polyval(para,xx);
        plot(xx,Y,'k-','linewidth',2);
        
        xlabel('Local Clustering Index');
        ylabel(sprintf('|\\Delta\\mu (\\circ)|'));
        title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))]);
        SetFigure(20);
        
        % ----  Plot 2D:  dPSE v.s. (dPrime & Clustering) -----
        
        figure(90 + fig_number); set(gcf,'Name',[area_list{a} ',' stim_type_list{k} ', dPrime & Clustering']); clf
        fig_number = fig_number + 1;
        
        scatter(LocCluster_S,dPrime_S,abs(dPSE(AMP_20_S,k)).^2*30+100,'ko','markerfacecolor','k');
        hold on;  scatter(LocCluster_NS,dPrime_NS,abs(dPSE(AMP_20_NS,k)).^2*30+100,'ko','Linewidth',2);
        set(gca,'yscale','log');
        
        set(findall(gca,'color','k'),'linewidth',2);
        set(gca,'ticklength',[0.02 0]);        
        
        xlabel('Local Clustering Index');
        ylabel('Discriminability (d'')');

        SetFigure(20);
         
        % -------- Data output
        dPSEvsTuning{a,k} = {[dPSE(AMP_20_S,k) dPrime_S LocCluster_S], [dPSE(AMP_20_NS,k) dPrime_NS LocCluster_NS]};
%}
       
    end
end

keyboard









