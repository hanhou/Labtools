%%
visOnCent = -1.192 * 1000;
visOnSD = 0.031968234 * 1000;

t = a(:,1) * 1000;

set(figure(),'color','w');
set(0,'defaultaxesxtickmode','auto');
plot(t,smooth(a(:,2),5),'k',t,smooth(a(:,3),5),'k--','LineWidth',2.5);
hold on;
xlim([min(t) max(t)]);
xlabel('Time to saccade onset (ms)');
ylabel('Firing rate (Hz)');
lims = axis;
line([0 0],lims(3:4),'color','k','linewidth',2);
line([150 150],lims(3:4),'color','k','linestyle','--');

% Substract background firing
% subDur = 1:find(t > visOnCent,1);
% a(:,2:end) = a(:,2:end) - repmat(mean(a(subDur,2:end),1),size(a,1),1);

% Instead of using transparency, we could create the patch opaque and
% behind the plot line, by specifying a very small negative z
p = patch([visOnCent-2*visOnSD visOnCent+2*visOnSD visOnCent+2*visOnSD visOnCent-2*visOnSD],...
    [lims(3) lims(3) lims(4) lims(4)],-eps*ones(1,4),[1 .9 .9],'EdgeColor','none');
%set(p,'FaceAlpha',0.5);
p = patch([visOnCent-2*visOnSD visOnCent+2*visOnSD visOnCent+2*visOnSD visOnCent-2*visOnSD]+1000,...
    [lims(3) lims(3) lims(4) lims(4)],-eps*ones(1,4),[.9 .9 .9],'EdgeColor','none');

%line([visOnCent-2*visOnSD visOnCent-2*visOnSD],lims(3:4),'color','k','linestyle','--');
%line([visOnCent+2*visOnSD visOnCent+2*visOnSD],lims(3:4),'color','k','linestyle','--');
% line([visOnCent-2*visOnSD visOnCent-2*visOnSD]+1,lims(3:4),'color','k','linestyle','--');
% line([visOnCent+2*visOnSD visOnCent+2*visOnSD]+1,lims(3:4),'color','k','linestyle','--');

axis(lims);
set(findall(gcf,'FontSize',10),'FontSize',15);

%%
colors = {'r','b','g','c','m'};
n_headings = size(a,2)-3;
set(figure(),'color','w');

for i = 1 : size(a,2)-3
    
    if i < (n_headings+1)/2     % Preferred headings
        colorSpec = ['-' colors{i}];
    elseif i == (n_headings+1)/2
        colorSpec = '-k';  % Zero heading
    else
        colorSpec = ['--' colors{n_headings-i+1}];    % Null headings
    end
    
    plot(t,smooth(a(:,i+3),20),colorSpec,'LineWidth',2);
    hold on;
    
end

xlim([min(t) max(t)]);
lims = axis;
line([0 0],lims(3:4),'color','k','linewidth',2);
line([150 150],lims(3:4),'color','k','linestyle','--');
xlabel('Time to saccade onset (ms)');
ylabel('Firing rate (Hz)');


% Instead of using transparency, we could create the patch opaque and
% behind the plot line, by specifying a very small negative z
p = patch([visOnCent-2*visOnSD visOnCent+2*visOnSD visOnCent+2*visOnSD visOnCent-2*visOnSD],...
    [lims(3) lims(3) lims(4) lims(4)],-eps*ones(1,4),[1 .9 .9],'EdgeColor','none');
%set(p,'FaceAlpha',0.5);
p = patch([visOnCent-2*visOnSD visOnCent+2*visOnSD visOnCent+2*visOnSD visOnCent-2*visOnSD]+1000,...
    [lims(3) lims(3) lims(4) lims(4)],-eps*ones(1,4),[.9 .9 .9],'EdgeColor','none');

%line([visOnCent-2*visOnSD visOnCent-2*visOnSD],lims(3:4),'color','k','linestyle','--');
%line([visOnCent+2*visOnSD visOnCent+2*visOnSD],lims(3:4),'color','k','linestyle','--');
% line([visOnCent-2*visOnSD visOnCent-2*visOnSD]+1,lims(3:4),'color','k','linestyle','--');
% line([visOnCent+2*visOnSD visOnCent+2*visOnSD]+1,lims(3:4),'color','k','linestyle','--');

axis(lims);
set(findall(gcf,'FontSize',10),'FontSize',15);


%%  Draw time course of CP and neurothreshold
psyThres = 2.7;

offset = visOnCent + 500;
set(figure(),'color','w');
[ax,h1,h2] = plotyy(b(:,1)+offset,b(:,2),b(:,1)+offset,b(:,3));
set(h1,'linewidth',2,'marker','o','markersize',10)
set(h2,'linewidth',2,'marker','o','markersize',10)
set(ax,'xlim',[b(1,1)+offset,b(end,1)+offset]);
set(get(ax(1),'Ylabel'),'String','Neuronal Threshold')
set(get(ax(2),'Ylabel'),'String','CP')

axes(ax(1));
set(gca,'yscale','log','ytickmode','auto');
lims = axis; hold on;
plot([0 0],[1 lims(4)],'k'); plot([lims(1) lims(2)],[psyThres psyThres],'--b');
ylim([2 lims(4)]);

axes(ax(2)); hold on;
set(plot(b(b(:,4)<0.05,1)+offset,b(b(:,4)<0.05,3),'o'),'markerfacecolor',get(h2,'color'),'markersize',10);   % Sig. CP
plot([lims(1) lims(2)],[0.5 0.5],'color',get(h2,'color'),'linestyle','--');
ylim([0.25 1]);
set(gca,'ytick',0.3:0.2:1);

set(findall(gcf,'FontSize',10),'FontSize',15);
xlabel('End of 500 ms time window (ms)');


%% Draw preferred directions
% Pref, DDI
a=[
    -174.21	0.787895
    62.1454	0.833732
    -79.4374	0.810634
    52.5521	0.891957
    -162.346	0.634979
    -151.382	0.774943
    81.8907	0.772405
    -139.346	0.808338
    -136.831	0.647282
    -136.363	0.626902
    -124.163	0.753386
    -104.249	0.700876
    -102.32	0.644762
    -96.9907	0.813763
    -93.6152	0.666371
    52.7113	0.877954
    78.8816	0.82837
    -83.7452	0.616777
    -69.7544	0.79881
    -67.613	0.774092
    -64.6874	0.74836
    -62.1465	0.787056
    -59.1189	0.736162
    -57.1985	0.531433
    46.8947	0.858053
    -55.8165	0.680515
    -50.7561	0.595275
    -43.2895	0.716236
    -34.1744	0.744204
    58.7105	0.816783
    59.5552	0.689749
    60.9656	0.723873
    -14.8404	0.728029
    -116.866	0.760786
    -23.227	0.720047
    -0.672953	0.762829
    12.8662	0.738146
    14.8967	0.65012
    22.658	0.616267
    29.3044	0.572263
    58.4656	0.788894
    -109.181	0.727989
    -118.263	0.822508
    71.5309	0.729416
    74.1824	0.796278
    70.5564	0.777844
    68.1016	0.790925
    87.587	0.754053
    96.3656	0.739684
    -118.058	0.706158
    62.3848	0.623138
    36.91	0.681146
    76.0593	0.762806
    62.7277	0.588263
    64.1009	0.687345
    -89.8141	0.685812
    -88.1316	0.695178
    37.6128	0.74175
    -59.4687	0.625262
    -98.4971	0.698903
    64.919	0.756944
    63.8811	0.692941
    -75.1856	0.699403
    -113.123	0.700464
    50.0637	0.711759
    47.0641	0.737244
    -32.3518	0.70059
    -90.9116	0.747768
    66.3604	0.582465
    44.7549	0.709515
    66.4095	0.722393
    34.5919	0.680917
    62.2171	0.676858
    -66.9159	0.728117
    -127.792	0.558973
    -49.972	0.811468
    -83.1822	0.731413
    62.7326	0.64371
    -103.052	0.692966
    -58.4016	0.712343
    -81.5543	0.701458
    61.1753	0.71318
    -83.4571	0.662489
    -115.82	0.742517
    -111.381	0.649143
    158.35	0.701732
    -122.867	0.630301
    67.6097	0.695042
    72.1934	0.550068
    -56.9986	0.625367
    -84.5979	0.734228
    15.7033	0.656405
    -120.988	0.629268
    -96.4109	0.64257
    43.0666	0.613104
    -107.552	0.598549
    -121.659	0.725204
    -122.902	0.64673
    13.5678	0.624464
    -90.9057	0.739774
    71.0705	0.677155
    74.3617	0.731458
    -58.8843	0.578604
    89.2469	0.688369
    85.7931	0.800397
    91.3107	0.771078
    101.375	0.834032
    126.084	0.520358
    -116.757	0.64542
    87.4187	0.556266
    -7.83134	0.472835
    -6.77033	0.574007
    -30.4597	0.547992
    -19.9005	0.679
    -113.998	0.468592
    -98.7648	0.558872
    -84.4403	0.520898
    155.3	0.706636
    174.956	0.455159
    28.553	0.48282
    -136.458	0.518638
    -77.5292	0.483265
    -46.0274	0.499357
    
    
    ];

set(0,'defaultaxesxtickmode','auto');
a = munique(a);
[~,ind]=sort(a(:,2),'descend');
a = a(ind,:);


% From heading to polar
a(:,1) = mod(90-a(:,1),360);


figure(11);roseInHeading(a(:,1)/180*pi);

figure(12);
set(12,'color','w');

% figure; hist(a(:,1),20)
for i = 1: size(a,1)
    polarwitherrorbar([ 0 a(i,1)/180*pi],[0 a(i,2)],[0 0],'b');
    hold on;
end
set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'fontsize',10),'fontsize',14);

%% d90 v.s. d270
figure;
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

% means = mean(abs(b(:,5:6)));
% errs = std(abs(b(:,5:6)))/sqrt(size(b,1));
%
% plot(means(1),means(2),'r*','markersize',10);
% plot([means(1)-errs(1) means(1)+errs(1)],[means(2) means(2)],'r-','linewidth',2);
% plot([means(1) means(1)], [means(2)-errs(2) means(2)+errs(2)],'r-','linewidth',2);


axis square;
[h p ] = ttest(abs(a(:,5))-abs(a(:,6)),0,0.05,'right')
% [h p ] = ttest(abs(b(:,5))-abs(b(:,6)),0,0.05)%,'right')

title(['p=' num2str(p)]);
figure();
plot(a(:,1),abs(a(:,5))./abs(a(:,6)),'o'); hold on;
xlim([0 360]);
plot([0 360],[0 0],'k-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitting tuning curve

%% Read xls data
[~,~,raw]=xlsread('results.xlsx',5);
%%
clear tuning16;
clear tuning8;
for i = 1:length(raw)
    if isnan(raw{i,:}) continue; end
    temp = str2num(raw{i,:});
    eval(['tuning' num2str(length(temp)) '(' num2str(i) ',:)= temp;']);
    i
end
tuning16 = munique(tuning16);
tuning16 = tuning16(sum(tuning16,2)>0,:);

tuning8 = munique(tuning8);
tuning8 = tuning8(sum(tuning8,2)>0,:);


%% Fitting

clear prefs sigma;
t = tuning8;
xx = (0:360/size(t,2):358)'*pi/180;
t = munique(t);
set(figure(),'color','w');

for i =  1:length(t)
    
    tt=t(i,:);
    subplot_tight(floor(sqrt(length(t))),ceil(length(t)/floor(sqrt(length(t)))),i,[0.01 0.01]);
    
    tried = 0;
    cf_.sigma = 1000;
    while tried < 10 && cf_.sigma > 2 * pi
        fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Algorithm','T', 'startpoint',[30,10,0.48,0.6],...
            'StartPoint',[10 10 3.1400000000000001 1.5 ],'Lower',[0 0 0 10*pi/180],'Upper',[Inf Inf   Inf Inf],'MaxFunEvals',5000);
        ft_ = fittype('a*exp(-2*(1-cos(x-pref))/sigma^2)+b',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a', 'b', 'pref', 'sigma'});
        [cf_, gof]= fit(xx,tt',ft_,fo_);
        
        tried = tried + 1
        confint(cf_)
        gof
    end
    
    plot(xx,tt); hold on;
    prefs(i) = cf_.pref;
    sigma(i) = cf_.sigma;
    
    if cf_.sigma < 2*pi
        plot(cf_);
        xlim([0 2*pi]); ylims = ylim();
        plot([pi/2 pi/2],[ylims(1) ylims(2)],'k--');
        plot([pi/2*3 pi/2*3],[ylims(1) ylims(2)],'k--');
    end
    
    ylims = ylim;
    
    text(0,ylims(2)*0.9,num2str(i));
    
    set(gca,'xtick',[],'ytick',[]);   xlabel(''); ylabel(''); legend off;
    drawnow;
    
end

set(findall(gcf,'color','b'),'linewidth',1.5,'marker','.');
set(findall(gcf,'color','r'),'linewidth',2);

% save('Tuning','prefs','sigma','tt');


%%
figure;
for i = 1: size(MU_prefs,2)
    polar([ 0 MU_prefs(i)],[0 1],'r');
    hold on;
end


for i = 1: size(SU_prefs,2)
    polar([ 0 SU_prefs(i)],[0 1]);
    hold on;
end


set(findall(gcf,'type','line'),'linewidth',2);
set(findall(gcf,'fontsize',10),'fontsize',14);


%%
clear
load tuning;

SU_prefs = SU_prefs*180/pi;
SU_sigma = SU_sigma*180/pi;

SU_prefs = -(SU_prefs-90);  % Transform to Gu's convention
SU_prefs = mod(SU_prefs + 180,360) - 180;  % To -180 ~ 180 range

MU_prefs = MU_prefs*180/pi;
MU_sigma = MU_sigma*180/pi;

MU_prefs = -(MU_prefs-90);  % Transform to Gu's convention
MU_prefs = mod(MU_prefs + 180,360) - 180;  % To -180 ~ 180 range


figure;
axes('position',[0.13 0.1 0.8 0.6]);

plot(MU_prefs,MU_sigma,'ro','markersize',10,'markerfacecolor','none'); hold on;
plot(SU_prefs,SU_sigma,'bo','markersize',10,'markerfacecolor','b');

set(gca,'xtickmode','auto')
axis([-180 180 0 270])
set(gca,'xtick',[-180 -90 0 90 180])
xlabel(sprintf('Preferred direction (\\circ)'));
ylabel(sprintf('Tuning width (\\circ)'));

set(gcf,'color','w');
set(findall(gcf,'fontsize',10),'fontsize',15)

axes('position',[0.13 0.7 0.8 0.2]);
box off;
hist([SU_prefs MU_prefs],30);
xlim([-180 180]);
set(gca,'xtick',[],'ytick',[]);
set(findobj(gca,'type','patch'),'edgecolor','k','facecolor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw spike waveform
%% Read xls data
[~,~,raw]=xlsread('results.xlsx',5);
clear spikeWF;
%%
figure; c=1; dt = 1/25; % ms
for i = 1:length(raw)
    if isnan(raw{i,:}) continue; end
    
    figure(1);
    spikeWF{c} = str2num(raw{i,:});
    maxPos = find(spikeWF{c}==1,1);
    trough1 = find(spikeWF{c}(1:maxPos)==min(spikeWF{c}(1:maxPos)));
    trough2 = find(spikeWF{c}(maxPos:end)==min(spikeWF{c}(maxPos:end))) + maxPos;
    up(c) = maxPos - trough1;
    down(c) = trough2 - maxPos;
    halfWid(c) = trough2 - trough1;
    % halfWid(c) = abs(find(spikeWF{c}==1,1)-find(spikeWF{c}==0));
    if halfWid(c)*dt<0.66 color = 'r'; else color='b'; end
    plot((0:dt:(length(spikeWF{c})-1)*dt)-maxPos*dt,(spikeWF{c}),color); hold on;
    figure(2);
    plot(spikeWF{c}(1:end-1),diff(spikeWF{c})); hold on;
    
    c=c+1;
end

% figure;
% plot(up,down,'o')
figure;
hist(halfWid*dt,30);

xlim([0.3 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CP
a=[
    6.9177	0.413008	0.071
    10.3393	0.432035	0.347
    43.3631	0.470588	0.769
    64.1971	0.496528	0.957
    201.808	0.579861	0.371
    45.9165	0.59075 	0.127
    8.37739	0.617949	0.002
    9.85292	0.681586	0.048
    ];

figure();
semilogx(a(:,1),a(:,2),'ko','markersize',10);
hold on;
semilogx(a(a(:,3)<0.05,1),a(a(:,3)<0.05,2),'ko','markerfacecolor','k','markersize',10);
set(line([1 1000],[0.5 0.5]),'linestyle',':','color','k');
[r,p]=corrcoef(log(a(:,1)),a(:,2))

[para,S]=polyfit(log(a(:,1)),a(:,2),1);
xx=3:400;
[Y,delta]=polyconf(para,log(xx),S);
plot(xx,Y,'r');
plot(xx,Y+delta,'r--',xx,Y-delta,'r--');

text(50,0.7,sprintf('\\itr \\rm= %3.2g\n\\itp \\rm= %3.2g',r(2),p(2)),'fontsize',20,'color','r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot thresholds
a=[1.97156	2.90268E-06	45.9165
    1.11736	3.26117E-06	43.3631
    1.11736	1.42792E-05	9.85292
    1.32516	2.08098E-06	64.1971
    1.32516	8.69461E-07	201.808
    1.93601	-3.20256E-07	10.3393
    2.09467	1.02902E-06	12.0297
    2.98288	8.26E-06	15.7538
    3.07855	-8.23E-06	6.9177
    2.07033	5.91E-07	20.4029
    2.07033	-1.13E-05	9.78095
    3.64472	0.000425546	4.6476
    2.85327	1.52355E-05	8.37739
    ];

figure();
loglog(a(:,1),a(:,3),'ko','markerfacecolor','k','markersize',7); hold on;
axis([.5 300 .5 300]);
set(line([.5 300],[.5 300]),'color','k');
axis square;
set(gca,'xtick',[1 10 100],'xticklabel',[1 10 100],'yticklabel',[1 10 100],'ytick',[1 10 100])
xlabel(sprintf('Psychophysical threshold (\\circ)'));
ylabel(sprintf('Neuronal threshold (\\circ)'));
SetFigure(15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw NC-SC relationship
a=[-0.1289	0.1438	-0.1054	0.6082
    -0.0435404	0.584593	-0.413975	0.110922
    -0.0406669	0.720224	0.00383938	0.988741
    0.064678	0.568684	-0.24261	0.365278
    0.0706382	0.53352	0.790452	0.000267818
    0.0889869	0.43248	-0.0755681	7.81E-01
    0.1059	0.2302	-0.1502	0.4639
    0.12407	1.18E-01	-0.320842	0.225658
    0.239764	0.0321845	0.217331	0.41878
    0.2418	0.0056	0.3026	0.1329
    0.242093	0.0304983	0.736632	0.00113524
    0.267419	0.0164816	0.904142	1.54E-06
    0.334754	0.00225354	0.856757	2.25E-05
    0.337903	0.00217356	0.394043	0.130993
    0.355192	0.00122479	0.525943	0.0363876
    0.534374	3.30E-07	0.568177	0.0216643
    ];

MST = 1;
MT = 2;
LIP = 3;

b=[LIP
    MST
    LIP
    LIP
    MST
    LIP
    LIP
    MST
    MST
    LIP
    MT
    MST
    MT
    MST
    MST
    LIP
    ];

colors=['b','r','g'];
figure;
for type = MST:LIP
    select = (b==type)&(a(:,2)<0.05);
    plot(a(select,3),a(select,1),[colors(type) 'o'],'markersize',10,'markerfacecolor',colors(type)); hold on;
    select = find((b==type)&(a(:,2)>=0.05));
    plot(a(select,3),a(select,1),[colors(type) 'o'],'markersize',10,'markerfacecolor','none','linewidth',2);
end

[r,p]=corrcoef(a(:,3),a(:,1))

[para,S]=polyfit(a(:,3),a(:,1),1);
xx = -1:0.1:1;
plot(xx,polyval(para,xx),'k--','linewidth',2);
[para,S]=polyfit(a(b==MST,3),a(b==MST,1),1);
plot(xx,polyval(para,xx),'b','linewidth',2);


xlabel(sprintf('\\bfr\\rm_{signal}'));
ylabel(sprintf('\\bfr\\rm_{noise}'));


lims = [-1 1 -.4 0.7];
axis(lims);
plot([0 0],[lims(3) lims(4)],'k--');
plot([lims(1) lims(2)],[0 0],'k--');
set(gca,'xtick',-1:0.5:1,'ytick',-0.4:0.2:0.6);

box off;
set(gca,'tickdir','o');
text(-0.9,0.4,sprintf('\\itr \\rm= %3.2g\n\\itp \\rm= %3.2g',r(2),p(2)),'fontsize',15,'color','k');
SetFigure(20);

%% Stimulus-dependent NC
a=[ -0.11175   -0.019656     0.55668    -0.23568    -0.34738      0.2521    -0.11474    -0.31492    0.046901    -0.39215    0.013544    0.069483   -0.029077     0.44053      0.3874    -0.54404
    0.1855     0.34483   0.0072901    -0.17206    -0.43301    -0.73321     0.46967     0.34977     0.52989     0.32353     0.68821     0.60985    -0.55574     0.88573     0.86958     0.64466
    0.22438   -0.060193    -0.57109    -0.58992     0.57156    -0.33473    -0.19701    0.091237    -0.58333    -0.13558      0.5905    -0.28456    -0.90322     0.41603    -0.11237    -0.78657
    0.32784     0.43333      0.1452     0.13687   -0.089897   -0.035742     0.40921    -0.21438     0.15904     0.10874     0.10588    0.090699    -0.22447     0.43902     0.54658   -0.090042
    0.40043     0.63116     0.70642     0.22286     0.55861     -0.3225     0.11392     0.31471     0.16841     0.37484     0.12308     0.81514     0.86275    -0.33386     0.38077     0.14032
    -0.44584     0.24461    -0.12656    -0.92758    -0.19057     0.41165    -0.33121    -0.44671       0.375     0.39193    0.023696    -0.31726   -0.077702     0.62103    -0.91962    0.036091
    -0.55574     0.77263     0.31766     0.15569      0.3439      0.3505    -0.41933     0.40967     0.80851     0.66169     0.86833     0.11689    -0.87411     0.77682    -0.79305     0.82806
    -0.61237 9.0649e-018     0.64624    -0.34719    -0.49383     0.27763     0.74227     0.92654     0.43298     0.93236     0.86483     0.91119    -0.11116     0.81692    -0.49828     0.24817
    0.68372     0.90434     0.33864     0.72488     0.17246     0.28653      0.2068    -0.75257     0.37326     0.21315       0.792     0.90758     0.87489     0.96483     0.77191     0.18171
    0.7441    0.099602    -0.30747    -0.71611      0.9268     0.10578    0.067416     -0.1743     0.10662     0.63595    -0.19361    -0.33602     0.50332     0.27609     0.92108     0.46771
    0.74592    -0.75366     0.64607     0.60075    0.094749     0.49031   -0.011156      0.3921     0.29644    -0.92859     0.71592     0.28448    -0.42848     0.24408     0.25138     0.60967
    -0.83366      0.6784    0.031793    -0.42535    -0.45691     0.58333     0.42601    -0.82391    0.084476     0.13558     0.50262    -0.73793     0.61201     0.38966    -0.13401    -0.38822
    0.90752    -0.53657    -0.21355     0.57498     0.33643     0.93763     0.16483     0.30189     0.41046     0.63828     0.65943     0.26558    0.041735     0.92997     0.17434     0.56976
    ];
PD1 = [31.5344
    221.093
    57.8543
    323.916
    8.83614
    75.1033
    149.119
    104.84
    152.258
    229.589
    93.3202
    57.8543
    104.84
    ];

PD2 = [173.745
    265.491
    265.491
    147.198
    29.6323
    70.355
    145.816
    149.119
    156.318
    150.461
    254.761
    221.093
    145.816
    ];

MST = 1;
MT = 2;
LIP = 3;

area = [MST
    LIP
    LIP
    MST
    MT
    MST
    MST
    MST
    LIP
    MST
    MT
    LIP
    MST
    ];

theta = 0:22.5:350;
inside = []; outside = [];

for i = 1:length(area)
    
    ind1 = min(find(PD1(i)<theta,1),find(PD2(i)<theta,1));
    ind2 = max(find(PD1(i)<theta,1),find(PD2(i)<theta,1));
    
    if abs(PD1(i)-PD2(i)) < 180
        inside = [inside a(i,ind1:ind2-1)];
        outside = [outside a(i,1:ind1-1) a(i,ind2:end)];
    else
        inside = [inside a(i,1:ind1-1) a(i,ind2:end)];
        outside = [outside a(i,ind1:ind2-1)];
    end
end

[h,p] = ttest2(inside,outside)
mean(inside)
mean(outside)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Microstimulation Pooled Data

clear;
MST = 1;
MT = 2;

ustim = [
    20	200	3	10	0.871	2.044	0.84	1.914	0.96 	0.82 	0.031	0.936399217	1	1
    100	200	3	10	0.7	2.468	0.567	2.404	0.82 	0.93 	0.133	0.974068071	1	1
    200	200	3	10	0.711	2.041	0.229	2.834	0.20 	0.28 	0.482	1.388535032	1	1
    20	200	2	10	-1.761	2.368	-2.278	3.014	0.96 	0.41 	0.517	1.272804054	1	0
    20	200	2	8	0.461	1.747	0.34	1.542	0.86 	0.70 	-0.121	0.882655982	1	0
    100	200	2	10	-0.03	1.186	0.707	1.816	0.08 	0.17 	0.737	1.531197302	1	0
    20	200	2	10	-0.059	1.011	-0.555	1.236	0.12 	0.50 	-0.496	1.222551929	0.5	0
    200	200	2	10	-0.55	1.084	-0.437	1.187	0.59 	0.76 	0.113	1.09501845	0.5	0
    50	200	2	10	0.058	2.234	-0.266	2.012	0.45 	0.73 	0.324	0.900626679	-1	0
    200	200	2	10	0.281	1.396	-1.078	2.251	0.00 	0.12 	1.359	1.612464183	-1	0
    50	200	2	9	-0.197	1.28	-0.221	2.097	0.84 	0.13 	-0.024	1.63828125	0.5	0
    200	200	2	10	0.688	1.222	-0.011	3.782	0.01 	0.00 	-0.699	3.09492635	0.5	0
    50	200	2	10	-0.654	1.652	-1.986	2.431	0.08 	0.20 	1.332	1.471549637	1	0
    20	200	2	10	-0.433	1.347	-1.336	2.781	0.48 	0.03 	0.903	2.064587973	1	0
    100	200	2	10	0.101	1.535	-1.021	1.419	0.00 	0.79 	1.122	0.924429967	1	0
    20	200	2	10	0.071	1.371	0.884	2.563	0.18 	0.05 	-0.813	1.869438366	1	0
    20	200	2	10	-0.817	1.754	-1.019	2.402	0.85 	0.29 	-0.202	1.369441277	1	1
    100	200	2	10	0.293	1.47	0.636	0.883	0.11 	0.13 	0.343	0.600680272	1	1
    100	200	2	10	1.952	1.934	3.268	3.673	0.64 	0.04 	1.316	1.899172699	0	0
    20	200	2	10	1.729	2.112	-0.48	2.212	0.00 	0.88 	2.209	1.047348485	1	0
    100	200	2	10	2.446	1.717	2.833	2.551	0.31 	0.19 	-0.387	1.485730926	1	0
    
    ];

tuning = [
    50.7791	0.679037	0.266983	2E-08	2.17195	1.02096	30.1622	76.3134
    50.7791	0.679037	0.266983	2E-08	2.17195	1.02096	30.1622	76.3134
    50.7791	0.679037	0.266983	2E-08	2.17195	1.02096	30.1622	76.3134
    56.6218	0.795531	0.743287	9E-42	4.73076	0.548376	48.336	102.202
    -106.969	0.722219	0.757427	3E-10	-3.07246	-5.36752	-109.785	121.656
    -106.969	0.722219	0.757427	3E-10	-3.07246	-5.36752	-109.785	121.656
    -70.2538	0.484358	0.603798	2E-04	-1.92213	0.0293397	-47.5109	104.056
    -70.2538	0.484358	0.603798	2E-04	-1.92213	0.0293397	-47.5109	104.056
    26.0371	0.487313	0.728425	9E-03	0.345092	0.330644	21.2946	178.858
    26.0371	0.487313	0.728425	9E-03	0.345092	0.330644	21.2946	178.858
    -72.2577	0.722165	0.481909	2E-06	-2.68169	-7.41287	-71.3109	150.512
    -72.2577	0.722165	0.481909	2E-06	-2.68169	-7.41287	-71.3109	150.512
    106.611	0.772796	0.650368	3E-12	3.7943	5.64588	104.356	202.82
    106.611	0.772796	0.650368	3E-12	3.7943	5.64588	104.356	202.82
    106.611	0.772796	0.650368	3E-12	3.7943	5.64588	104.356	202.82
    97.3075	0.58952	0.65222	6E-09	4.89053	0.951365	90	139.263
    -89.556	0.541205	0.319884	3E-02	-2.51895	-1.24304	-86.5258	169.122
    -89.556	0.541205	0.319884	3E-02	-2.51895	-1.24304	-86.5258	169.122
    -69.3562	0.36755	0.213621	8E-01	-0.314286	0.165432	-45.0032	11.6246
    30.8287	0.826801	0.68124	3E-19	5.22725	0.671751	27.0359	59.1604
    30.8287	0.826801	0.68124	3E-19	5.22725	0.671751	27.0359	59.1604
    
    ];

%%
load ustim_MST;
% load ustim_VIP;

amp = ustim(:,1);
reps = ustim(:,4);
NoSigma = ustim(:,3);
pPSE = ustim(:,9);
pSIG = ustim(:,10);
dPSE = ustim(:,11);
dSIG = ustim(:,12);
tuningOK = ustim(:,13);
electrodeFailure = ustim(:,14);


%% PSE shift


figure(91); clf

%%%%%% Mask here...

% 20140512 MST_20
% mask_PSE = (0|reps >= 5) & (0|tuningOK > 0) & (0|amp <= 20) & ~(electrodeFailure);

% 20140512 MST_all
mask_PSE = (0|reps >= 5) & (0|tuningOK > 0.5) & (1|amp <= 50) & ~(electrodeFailure);

S_PSE = (pPSE < 0.05) & mask_PSE;
NS_PSE= (pPSE >= 0.05) & mask_PSE;

bin = 0.5;
xcenters = round(min(dPSE(S_PSE|NS_PSE))/2)*2-bin/2:bin:round(max(dPSE(S_PSE|NS_PSE))/2)*2+0.5;  % Ensure that 0 is one of the borders.

%%%%%%


histS = hist(dPSE(S_PSE),xcenters);
histNS = hist(dPSE(NS_PSE),xcenters);
hbars = bar(xcenters,[histS' histNS'],1,'stacked','k','LineWidth',2);
set(hbars(2),'FaceColor','none');

all_PSE = dPSE(S_PSE|NS_PSE);
meanPSE = mean(all_PSE);
maxY = max([histS+histNS])*1.2;
axis([min(xcenters)-1 max(xcenters)+1 0 maxY]);
set(gca,'xtick',[-100:2:100],'xMinorTick','on');
xlim([min(-4,min(all_PSE))-0.5 max(10,max(all_PSE))+0.5]);

% ttest
[~,p] = ttest(all_PSE,0);

hold on; plot([meanPSE meanPSE],[maxY/1 maxY/1.2],'k'); fill([meanPSE meanPSE+0.2 meanPSE-0.2],[maxY/1.2 maxY/1.2+0.5 maxY/1.2+0.5],'k');
text(meanPSE+1, maxY/1.2+0.5, sprintf('Mean \\DeltaPSE = %2.2f^\\circ\n\\itp \\rm= %.2g\nN = %g',meanPSE,p,sum(S_PSE)+sum(NS_PSE)));

plot([0 0],[0 maxY],'k:' );

xlabel(sprintf('Induced PSE shift (\\circ)'));
ylabel('Number of cases');
legend p<0.05 p>0.05
SetFigure()

%% 2-D Plot with different conditions
figure(92); clf;
result_200 = [];
result_50_100 = [];
result_20_20 = [];
result_20_45 = [];

% all conditions
mask_2D = (1|reps >= 5) & (1|tuningOK > 0) & (1|amp <= 50) & ~(electrodeFailure);


for i = find(mask_2D)'
    
    if amp(i) <= 50
        %         if ~(tuningOK(i)==1); continue; end;
        markerSize = 10;
    elseif  amp(i) <= 100
        markerSize = 20;
        result_50_100 = [result_50_100; dPSE(i) dSIG(i)];
    elseif amp(i) >= 200
        markerSize = 40;
        result_200 = [result_200; dPSE(i) dSIG(i)];
    end
    
    if pPSE(i) < 0.05 && pSIG(i) < 0.05
        markerCol = 'm';
        %         markerFaceCol = 'm';
    elseif pPSE(i) <= 0.05 && pSIG(i) >= 0.05
        markerCol = 'r';
        %         markerFaceCol = 'r';
    elseif pPSE(i) >= 0.05 && pSIG(i) < 0.05
        markerCol = 'b';
        %         markerFaceCol = 'b';
    else
        markerCol = 'k';
        %         markerFaceCol = 'none';
    end
    
    markerFaceCol = 'none';
    
    if NoSigma(i) == 2.0 && amp(i) == 20
        %         marker = '^';
        marker = 'o';
        result_20_20 = [result_20_20; dPSE(i) dSIG(i)];
    elseif amp(i) == 20
        result_20_45 = [result_20_45; dPSE(i) dSIG(i)];
        marker = 'o';
    else
        marker = 'o';
    end
    
    plot(dPSE(i),dSIG(i),[markerCol marker],'markersize',markerSize,'linewidth',2,'markerFaceCol',markerFaceCol);
    set(gca,'yScale','log');
    hold on;
    
end

text(meanPSE+1, maxY/1.2+0.5, sprintf('\\itN\\rm = %g',sum(mask_2D)));


lims = axis;
plot([lims(1)*1.1 lims(2)*1.1],[1 1],'k:');
plot([0 0],[lims(3)/1.2 lims(4)*1.2],'k:');
axis([lims(1:2)*1.1 lims(3)/1.2 lims(4)*1.2]);
SetFigure(20);
xlabel(sprintf('\\DeltaPSE (\\circ)'),'color','r');
ylabel(sprintf('\\sigma_{ratio}'),'color','b');

% axis([-6 15 0.5 6]);

%% dPSE v.s. Discrimi
figure(93); clf;

% 20140512 MST_100uA
mask_Discri = (0|tuningOK > 0) & (0|amp <= 100);

% 20140512 VIP_20uA, d'0degree
% mask_Discri = (0|tuningOK > 0.5) & (0|amp <= 50);

AMP_20_S = mask_Discri & (pPSE<0.05);
AMP_20_NS = mask_Discri & (pPSE>=0.05);
% AMP_20_S = AMP_20_S(1:size(tuning,1));
% AMP_20_NS = AMP_20_NS(1:size(tuning,1));

%%{
% HTI and dPSE (linear)

% DDI
dis_S = abs(tuning(AMP_20_S,2));
dis_NS = abs(tuning(AMP_20_NS,2));

% Deviation of pref. heading from +/- 90
% dis_S = abs((abs(tuning(AMP_20_S,1))-90));
% dis_NS = abs((abs(tuning(AMP_20_NS,1))-90));


plot(dis_S,abs(dPSE(AMP_20_S)),'ro','markerfacecolor','r');
hold on;  plot(dis_NS,abs(dPSE(AMP_20_NS)),'ro');
set(findall(gca,'color','r'),'markersize',10);
set(findall(gca,'color','r'),'linewidth',2);
set(gca,'ticklength',[0.02 0]);

% fitting
allX = ([dis_S; dis_NS]);
allY = (abs([dPSE(AMP_20_S); dPSE(AMP_20_NS)]));

axis([min(allX)/1.5 max(allX)*1.5 min(allY)/1.5 max(allY)*1.5]);

[r,p]=corrcoef(allX,allY);
[para,S]=polyfit(allX,allY,1);
xx = min(allX)-0.5:0.1:max(allX)+0.5;
Y = polyval(para,xx);
plot(xx,Y,'r-','linewidth',2);

xlabel('DDI');
ylabel(sprintf('|Induced PSE shift (\\circ)|'));
title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))]);
SetFigure(20);

%}


%%{
% D-prime and dPSE (X log scale)
figure(94); clf;
%
% dis_S = (abs(tuning(AMP_20_S,6))+abs(tuning(AMP_20_S,5)))/2;
% dis_NS = (abs(tuning(AMP_20_NS,6))+abs(tuning(AMP_20_NS,5)))/2;
%
dis_S = abs(tuning(AMP_20_S,6));
dis_NS = abs(tuning(AMP_20_NS,6));

loglog(dis_S,abs(dPSE(AMP_20_S)),'ro','markerfacecolor','r');
hold on;  plot(dis_NS,abs(dPSE(AMP_20_NS)),'ro');
set(findall(gca,'color','r'),'markersize',10);
set(findall(gca,'color','r'),'linewidth',2);
set(gca,'ticklength',[0.02 0]);

% fitting
allX = ([dis_S; dis_NS]);
allY = (abs([dPSE(AMP_20_S); dPSE(AMP_20_NS)]));

axis([min(allX)/1.5 max(allX)*1.5 min(allY)/1.5 max(allY)*1.5]);

[r,p]=corrcoef(log(allX),log(allY));
[para,S]=polyfit(log(allX),log(allY),1);
xx = min(log(allX))-0.5:0.1:max(log(allX))+0.5;
Y = polyval(para,xx);
plot(exp(xx),exp(Y),'r-','linewidth',2);

xlabel('Discriminability (d'')');
ylabel(sprintf('|Induced PSE shift (\\circ)|'));
title(['r = ' num2str(r(2)) ', p = ' num2str(p(2))]);
set(gca,'xtick',[0.1 :0.1 1:1:10]);
SetFigure(20);

%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw tissue damage
x1 = [10198 10392 10586 10835 10929]' - 10571;
y1 = [0.2 0.5 1.5 0.5 0.2]';
x2 = [11295 11377 11410 11454 11520 11589 11639 11670 11720]' - 11620;
y2 = [0.2 0.3 0.7 1.0 1.1 1.1 0.4 0.2 0.2]';
x3 = [-530 -456 -422 -400 -300 -200 -100 0 100 200 226 250 265 289 299 372]';
y3 = [0.1 0.2 0.4 0.8 1 1.1 1.1 1.1 1.1 1.1 0.9 0.6 0.5 0.3 0.2 0.1]';

figure();
plot(x1,y1,'ok','markerfacecolor','k','markersize',10); hold on;
plot(x2,y2,'or','markerfacecolor','r','markersize',10);
xx = -500:500;

ylim([0 1.6])
plot([-400 400],[0.2 0.2],'k--');
plot([0 0],[0 1.6],'k--')

% Gaussian fitting
f1 = fittype('Gauss1');
p1 = fit(x1,y1-y1(1),f1);
plot(xx,p1(xx)+y1(1),'k','linewidth',2);

f2 = fittype('Gauss1');
p2 = fit(x2,y2-y2(1),f1);
plot(xx,p2(xx)+y2(1),'r','linewidth',2);

xlabel(sprintf('Distance to \\mu-stim site (\\mum)'));
ylabel('Impedance (AU)');



%%  Draw SLEEPING MONKEY  20131030

load sleeping.mat;
xcenters = 0:3:300;
sleepFR = hist(sleep,xcenters);
sleepFR = smooth(sleepFR / sum(sleepFR));
awakeFR = hist(awake,xcenters);
awakeFR = smooth(awakeFR / sum(awakeFR));
dreamFR = hist(dream,xcenters);
dreamFR = smooth(dreamFR / sum(dreamFR));
figure();
plot(xcenters,awakeFR,'r',xcenters,sleepFR,'k',xcenters,dreamFR,'b','linewidth',3);
hold on;
plot(median(sleep),0.04,'sk',median(awake),0.04,'sr',median(dream),0.04,'bs');
legend('Awake','Sleep', 'Dream');
xlabel('Spontaneous firing rate (Hz)');
ylabel('Probability');
title('Distribution of instantaneous MU firing rate during different brain states');

%% Draw sleeping NC
load sleepingNC.mat
cell1 = V20131029_SLEEPINGMONKEY_3_c115_state_trasition_Ch7;
cell2 = V20131029_SLEEPINGMONKEY_3_c115_state_trasition_Ch8;

periods{1} = [0 30; 80 112; 1218 1265; 1664 1713; 2003 2118; 2391 2474; 2540 2608 ];  % awake in s
periods{2} = [278 358; 651 773; 814 927;1292 1558;];  % sleep in s

colors = {'r','k'};

select = {[] []};

binWid = 0.5; % in s
t = 0:binWid:max(cell1.times);

cell1FR = hist(cell1.times,t)/binWid;
cell2FR = hist(cell2.times,t)/binWid;

% sorting to different period
for i = 1:2
    for j = 1:size(periods{i},1)
        select{i} = [select{i} find(t>periods{i}(j,1) & t<periods{i}(j,2))]; % awake
    end
end


for i = 1:2
    cell1FR_temp = cell1FR(select{i});
    cell2FR_temp = cell2FR(select{i});
    
    figure(2000+i); clf; plot(cell1FR_temp,cell2FR_temp,[colors{i} 'o'],'markerfacecolor',colors{i},'markersize',3);
    
    cell1FR_z_score{i} = (cell1FR_temp-mean(cell1FR_temp))/std(cell1FR_temp);
    cell2FR_z_score{i} = (cell2FR_temp-mean(cell2FR_temp))/std(cell2FR_temp);
    
    figure(2003+i); clf;  plot(cell1FR_z_score{i},cell2FR_z_score{i},[colors{i} 'o'],'markerfacecolor',colors{i},'markersize',3);
    [r,p] = corrcoef(cell1FR_z_score{i},cell2FR_z_score{i});
    
    para = polyfit(cell1FR_z_score{i},cell2FR_z_score{i},1);
    xmin = min([cell1FR_z_score{i},cell2FR_z_score{i}]);
    xmax = max([cell1FR_z_score{i},cell2FR_z_score{i}]);
    xlim([xmin, xmax]);
    ylim([xmin, xmax]);
    axis square;
    xlabel('Z-scored response (Cell 1)');
    ylabel('Z-scored response (Cell 2)');
    hold on;
    
    plot([xmin,xmax],polyval(para,[xmin,xmax]),[colors{i} '--'],'linewidth',3);
    %plot([xmin,xmax], [xmin, xmax],'k--');
    box on;
    
    title(['NC = ' num2str(r(2)) ', p = ' num2str(p(2))]);
    
end


%% 20140505  Draw HTI and DDI

data = [
    0.36755	0.213621	7.56E-01
    0.36755	0.213621	7.56E-01
    0.401186	0.265314	0.659617
    0.410916	0.705438	0.356219
    0.424534	0.0941637	0.970321
    0.434226	0.0884684	0.965637
    0.46641	0.635087	0.149092
    0.474502	0.0880181	0.721181
    0.5156	0.1313	0.2888
    0.515676	0.12173	0.210623
    0.517871	0.244384	0.601671
    0.528	0.1522	0.04943
    0.529563	0.566418	6.16E-06
    0.541205	0.319884	0.0268448
    0.541205	0.319884	0.0268448
    0.541205	0.319884	0.0268448
    0.546071	0.0211337	0.90607
    0.569237	0.645948	0.000970349
    0.5694	0.38	0.3368
    0.578254	0.190447	0.62452
    0.58952	0.65222	6.43E-09
    0.58952	0.65222	6.43E-09
    0.593718	0.0885783	0.750432
    0.59483	0.022991	0.929311
    0.594972	0.643188	0.0560567
    0.595713	0.580307	2.56E-05
    0.601077	0.354711	0.531745
    0.604542	0.757075	0.000168388
    0.60857	0.162925	0.38076
    0.615632	0.50436	0.407188
    0.615632	0.50436	0.407188
    0.617249	0.711029	0.000145974
    0.620321	0.694074	0.434044
    0.622305	0.0349241	0.0850509
    0.622305	0.0349241	0.0850509
    0.632841	0.182315	0.106788
    0.646629	0.302961	0.00182973
    0.651758	0.296821	0.437297
    0.660662	0.456279	0.00528251
    0.660821	0.210695	0.000264052
    0.660821	0.210695	0.000264052
    0.664836	0.477801	0.0136163
    0.671996	0.0685874	0.0475821
    0.672716	0.294947	7.36E-06
    0.675122	0.199345	0.297423
    0.675447	0.133153	0.292315
    0.6785	0.56	0.006079
    0.679055	0.196242	0.0490595
    0.679838	0.133833	0.0855625
    0.680552	0.737971	2.44E-09
    0.680552	0.737971	2.44E-09
    0.683143	0.648106	0.0163491
    0.688907	0.179893	1.11E-06
    0.698216	0.264324	0.0948593
    0.698723	0.0898791	0.346594
    0.699669	0.26348	0.000584924
    0.705383	0.744171	4.48E-09
    0.718171	0.570198	0.00277561
    0.719096	0.189573	0.0455281
    0.72221	0.341434	0.00175455
    0.723664	0.364597	0.0136791
    0.732012	0.0552027	0.00888653
    0.732012	0.0552027	0.00888653
    0.732012	0.0552027	0.00888653
    0.733091	0.601468	6.04E-06
    0.734754	0.231829	2.79E-07
    0.737614	0.0747185	0.00615718
    0.739532	0.0748713	0.105154
    0.740061	0.322524	0.0149806
    0.740194	0.166817	0.00250216
    0.745259	0.190558	0.00204098
    0.745259	0.190558	0.00204098
    0.745259	0.190558	0.00204098
    0.746362	0.123275	0.000857029
    0.752857	0.480665	6.14E-04
    0.757074	0.651415	5.64E-09
    0.76303	0.542112	2.22E-05
    0.766785	0.313565	2.06E-11
    0.769053	0.453743	4.24E-04
    0.774648	0.422414	0.00733996
    0.777269	0.318097	1.09E-08
    0.780203	0.208311	8.79E-06
    0.78162	0.594875	9.09E-07
    0.78162	0.594875	9.09E-07
    0.785803	0.248381	2.28E-07
    0.785901	0.719458	5.47E-06
    0.787232	0.714175	3.39E-05
    0.789513	0.605814	4.05E-07
    0.789513	0.605814	4.05E-07
    0.789513	0.605814	4.05E-07
    0.790071	0.280533	8.37E-07
    0.790071	0.280533	8.37E-07
    0.79139	0.409155	9.22E-05
    0.796319	0.67728	5.93E-08
    0.796319	0.67728	5.93E-08
    0.796319	0.67728	5.93E-08
    0.796319	0.67728	5.93E-08
    0.796319	0.67728	5.93E-08
    0.798671	0.107798	0.000281899
    0.798671	0.107798	0.000281899
    0.798713	0.487114	1.88E-05
    0.800588	0.115253	1.16E-05
    0.801273	0.461337	2.18E-10
    0.801273	0.461337	2.18E-10
    0.80483	0.632989	0
    0.80483	0.632989	0
    0.805729	0.542012	3.57E-09
    0.814266	0.346749	4.30E-10
    0.820928	0.625242	0
    0.821853	0.176034	1.84E-09
    0.823257	0.368382	1.97E-08
    0.825488	0.520493	1.26E-07
    0.826801	0.68124	2.75E-19
    0.826801	0.68124	2.75E-19
    0.826801	0.68124	2.75E-19
    0.827128	0.245489	3.04E-07
    0.830115	0.706078	2.39E-07
    0.830115	0.706078	2.39E-07
    0.830115	0.706078	2.39E-07
    0.830115	0.706078	2.39E-07
    0.830844	0.54078	1.75E-07
    0.830844	0.54078	1.75E-07
    0.832184	0.279439	2.32E-09
    0.833001	0.480977	6.92E-08
    0.833001	0.480977	6.92E-08
    0.833001	0.480977	6.92E-08
    0.833521	0.53861	1.31E-07
    0.833521	0.53861	1.31E-07
    0.834335	0.261255	2.02E-08
    0.837357	0.667732	3.44E-15
    0.837357	0.667732	3.44E-15
    0.837706	0.645092	3.22E-15
    0.839101	0.43629	0.00345124
    0.840259	0.346198	5.72E-07
    0.840259	0.346198	5.72E-07
    0.840983	0.571353	4.88E-15
    0.843483	0.512295	5.94E-07
    0.843691	0.495806	1.21E-09
    0.845493	0.74417	1.24E-07
    0.845493	0.74417	1.24E-07
    0.845493	0.74417	1.24E-07
    0.846154	0.509866	2.53E-05
    0.846842	0.78647	1.30E-14
    0.852808	0.606105	0
    0.852931	0.457084	6.07E-08
    0.854181	0.738912	1.61E-13
    0.854398	0.579856	1.32E-13
    0.854398	0.579856	1.32E-13
    0.854398	0.579856	1.32E-13
    0.854398	0.579856	1.32E-13
    0.854423	0.509363	3.97E-05
    0.854423	0.509363	3.97E-05
    0.855903	0.755424	1.11E-16
    0.857701	0.667242	1.51E-13
    0.858935	0.477505	3.78E-08
    0.859328	0.540286	1.11E-16
    0.859569	0.0960388	2.95E-09
    0.860296	0.548533	2.33E-15
    0.861236	0.348241	1.22E-09
    0.86238	0.70529	0
    0.864551	0.579017	5.22E-15
    0.865886	0.672367	9.95E-09
    0.867953	0.312605	3.66E-07
    0.868604	0.501578	0
    0.869683	0.664873	2.13E-10
    0.870636	0.236819	6.66E-16
    0.872252	0.734584	6.93E-08
    0.880209	0.402709	0
    0.881262	0.708563	0
    0.881324	0.581054	0
    0.882641	0.666835	8.66E-15
    0.886188	0.619019	0
    0.886188	0.619019	0
    0.886434	0.298877	7.02E-10
    0.888389	0.681339	0
    0.888389	0.681339	0
    0.888954	0.421133	4.60E-06
    0.89026	0.52893	0
    0.892368	0.377042	1.11E-11
    0.892431	0.738788	9.36E-12
    0.892431	0.738788	9.36E-12
    0.896402	0.618071	0
    0.898147	0.752359	3.24E-11
    0.899282	0.745891	1.33E-12
    0.899282	0.745891	1.33E-12
    0.899282	0.745891	1.33E-12
    0.899282	0.745891	1.33E-12
    0.900636	0.698979	0
    0.901425	0.362693	0
    0.90171	0.540624	0
    0.901937	0.918845	0
    0.910095	0.637755	3.15E-07
    0.912167	0.510502	1.60E-13
    0.912167	0.510502	1.60E-13
    0.912167	0.510502	1.60E-13
    0.912167	0.510502	1.60E-13
    0.912802	0.734927	9.99E-16
    0.913451	0.250924	0
    0.914915	0.68882	0
    0.914915	0.68882	0
    0.921719	0.697609	0
    0.932658	0.502984	0
    0.940557	0.565476	0
    0.946619	0.599037	0
    0.946619	0.599037	0
    
    ];

data = munique(data);
figure();
scatter(data(data(:,3)<0.05,1),data(data(:,3)<0.05,2),7 * -log(data(data(:,3)<0.05,3)+eps),'or');
hold on;
plot(data(data(:,3)>=0.05,1), data(data(:,3)>=0.05,2),'o');
axis square;
axis([0 1 0 1]);
xlabel('DDI');
ylabel('HTI');
SetFigure(15);




%% Transform function test in MOOGDots  20141031
t = 0:0.01:10;
model = diff(normcdf(t,0,1));
x = model + randn(size(model))*3;

% zs = 6.94525;
% ps = 46.08874;

zs = -20:.05:20;
ps = 0.1 * 1.1.^(1:100);

error_to_model = zeros(length(zs),length(ps));
error_to_data = zeros(length(zs),length(ps));
tic

for iz = 1:length(zs)
    z = zs(iz);
    
    for ip = 1:length(ps)
        p = ps(ip);
        
        y = x;
        for i = 2:length(model)
            y(i) = (1/(1+p)) * (-y(i-1)*(1-p)+x(i)*(1+z) + x(i-1)*(1-z));
        end
        
%         figure(91); clf;
%         plot(t,x,'ro-'); hold on
%         plot(t,y,'b.--');
%         text(1,-max(abs(x))*.8,sprintf('z = %g\np = %g',z,p),'fontsize',20);
%         ylim([-max(abs(x)) max(abs(x))]);
%         drawnow;

        error_to_model(iz,ip) = log(norm(y-model));
        error_to_data(iz,ip) = log(norm(y-x));
        
    end
end
toc

model_to_data = log(norm(model-x));

[zz,pp] = meshgrid(zs,log(ps));

figure(7); clf; 
surfc(zz,pp,error_to_model','Edgecolor','none');

xlabel('z'); ylabel('log (p)');  title('Error to model'); colorbar; view(2); axis tight; SetFigure();

figure(8); clf; surfc(zz,pp,error_to_data','Edgecolor','none');
xlabel('z'); ylabel('log (p)');  title('Error to data'); colorbar; view(2); axis tight;  SetFigure();
 


%% P-value, AUC, and d'-prime
du = -20:0.5:20;
sig = 0.1:0.2:10;

n = 20;
A = randn(n,length(du)*length(sig));
B = randn(n,length(du)*length(sig));
AUC = zeros(1,length(du)*length(sig));

i = 1;
for uu = du
    for ss = sig
        B(:,i) = B(:,i) * ss + uu; % Shift to (uu,ss)
        AUC(i) = abs(rocN(B(:,i),A(:,i))-0.5);
        i = i + 1;
    end
end

[~,ps] = ttest2(A,B);  % p-value 
sqrtlogps = sqrt(-log(ps));

dprime = abs((mean(A) - mean(B)))./sqrt((std(A).^2 + std(B).^2)/2); % d-prime


figure(99);  
clf;

plim = 5;

subplot(1,2,1);
hold on; 
plot(sqrtlogps(sqrtlogps<plim), dprime(sqrtlogps<plim),'r.');  axis tight;
plot([sqrt(-log(0.05)) sqrt(-log(0.05))],ylim,'k--'); text(sqrt(-log(0.05)),0.7,'0.05');
plot([sqrt(-log(1e-3)) sqrt(-log(1e-3))],ylim,'k--'); text(sqrt(-log(1e-3)),0.7,'0.001');
plot([sqrt(-log(1e-10)) sqrt(-log(1e-10))],ylim,'k--'); text(sqrt(-log(1e-10)),0.7,'1e-10');
plot([sqrt(-log(1e-20)) sqrt(-log(1e-20))],ylim,'k--'); text(sqrt(-log(1e-20)),0.2,'1e-20');
xlabel('Sqrt(-log(p))');
ylabel('dprime');

subplot(1,2,2);
hold on;
xx = sqrtlogps(sqrtlogps<plim)';
yy =  AUC(sqrtlogps<plim)';

plot(xx, yy,'r.');  axis tight;
plot([sqrt(-log(0.05)) sqrt(-log(0.05))],ylim,'k--'); text(sqrt(-log(0.05)),0.2,'0.05');
plot([sqrt(-log(1e-3)) sqrt(-log(1e-3))],ylim,'k--'); text(sqrt(-log(1e-3)),0.2,'0.001');
plot([sqrt(-log(1e-10)) sqrt(-log(1e-10))],ylim,'k--'); text(sqrt(-log(1e-10)),0.2,'1e-10');
plot([sqrt(-log(1e-20)) sqrt(-log(1e-20))],ylim,'k--'); text(sqrt(-log(1e-20)),0.2,'1e-20');

xlabel('Sqrt(-log(p))');
ylabel('AUC');


% Fitting

[r,p] = corr(xx,yy,'type','spearman');
[para,S]=polyfit(xx,yy,1);

xxx = min(xx):0.1:max(xx);
Y = polyval(para,xxx);
plot(xxx,Y,'k','linewidth',2);  hold on;

text(0,0,sprintf('r = %3.3f, p = %3.3f',r,p),'color','k');
text(0,0,sprintf('y = %3.3f * x + %3.3f',para(1),para(2)),'color','k');


SetFigure();


%% ===================================
% Probability of heading angle change (shift) in constant HD task.  HH20150115

simN = 20000;
xcenterN = 100;

unique_headingNs = 2:2:40;
repNs = 200;


means = zeros(length(unique_headingNs),length(repNs));
stds = zeros(length(unique_headingNs),length(repNs));
xcenters = zeros(length(unique_headingNs),length(repNs),xcenterN);
props = zeros(length(unique_headingNs),length(repNs),xcenterN);


for hh = 1:length(unique_headingNs)
    
    hh
    unique_headingN = unique_headingNs(hh);
    unique_heading = -fix(unique_headingN/2): fix(unique_headingN/2);
   
    if mod(unique_headingN,2) == 0
        unique_heading(unique_heading==0)=[];
    end
    
    for pp = 1:length(repNs)
        repN = repNs(pp);
        
        headings = zeros(unique_headingN,repN*simN);
        for i = 1:repN*simN
            headings(:,i) = unique_heading(randperm(unique_headingN));
        end
        
        headings = reshape(headings,unique_headingN*repN,[]);
        shift_prop = sum(sign(headings(1:end-1,:)) ~= sign(headings(2:end,:)))/(size(headings,1)-1);
        
        [n,x] = hist(shift_prop,xcenterN);
        
        xcenters(hh,pp,:) = x;
        props(hh,pp,:) = n / simN / (x(2)-x(1));
        
        means(hh,pp) = mean(shift_prop);
        stds(hh,pp) = std(shift_prop);
        
        figure(3);
        plot(x,n); drawnow;
    end
end

ys = squeeze(props(:,1,:))';
xs = squeeze(xcenters(:,1,:))';

%%
figure(999); clf; hold on;

ft = fittype('gauss1');

for i = 1:size(xs,2)
    [fitresult, gof] = fit( xs(:,i),ys(:,i), ft );
    
    fitresults{i} = fitresult;
    means_fit(i) = fitresult.b1;
    stds_fit(i) = fitresult.c1;
    
    plot(fitresult);
end

plot(xs,ys);

figure(1999); clf; hold on; 
errorbar(unique_headingNs,means,stds);
errorbar(unique_headingNs,means_fit,stds_fit,'r');










