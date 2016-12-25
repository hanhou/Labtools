function Group_Psychometric_Polo(XlsData)

num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;

%% Get data

n = size(num,1);

date = num(:,1);
sigmas = num(:,[4 6 8]);
us = num(:,[5 7 9]);
coh = num(:,3);
tf = num(:,2);
glass = num(:,11);

platform = num(:,12);
showPlatform = 1;
platforms = [301 109 103 -103];
platformColor = {'k','r','g','c'};


colors = {'b','r','g'};

sigma_pred = 1./sqrt(1./sigmas(:,1).^2 + 1./sigmas(:,2).^2);
u_pred = (sigmas(:,2).^2 .* us(:,1) + sigmas(:,1).^2 .* us(:,2))./(sigmas(:,1).^2 + sigmas(:,2).^2);

%% Sigma time course
set(figure(3),'position',[19 85 1388 535]); clf
for k = 1:3
    plot(sigmas(:,k),[colors{k} 'o'],'markerfacecol',colors{k},'markersize',10); hold on;
    plot(smooth(sigmas(:,k),30,'rloess'),['-' colors{k}],'linewid',2);
end
plot(sigma_pred,'go','markersize',10,'linewidt',1.5);
% plot(coh);
xlabel('Session');
ylabel('Threshold');

ylim([1 20]); xlim([0 n+2]);
set(gca,'YScale','log','Ytick',[1:10 20]); 

Annotation(date,tf,glass);

SetFigure();

%% u time course
set(figure(4),'position',[19 85 1388 535]); clf
plot([0 n],[0 0],'k--','linew',2); 
hold on

for k = 1:3
    plot(us(:,k),[colors{k} 'o'],'markerfacecol',colors{k},'markersize',10);
    plot(smooth(us(:,k),30,'rloess'),['-' colors{k}],'linewid',2);
end
% plot(u_pred,'go','markersize',13,'linewidt',1.5);
xlabel('Session');
ylabel('Bias');
ylim([-5 5]); xlim([0 n+2]);

Annotation(date,tf,glass);


SetFigure();

%% Prediction ratio time course
sigma_pred_ratio = sigmas(:,3)./sigma_pred;

set(figure(5),'position',[19 85 1388 535]); 
clf;
plot([0 n],[1 1],'k--','linew',2); hold on;

failureBegin = 72;

% Min threshold
% plot(min(sigmas(:,1:2),[],2)./sigma_pred,'v','color',[0.5 0.5 0.5],'markerfacecol',[0.5 0.5 0.5]);
plot(smooth(min(sigmas(:,1:2),[],2)./sigma_pred,20,'rloess'),['-.k'],'linewid',2);
% plot(max(sigmas(:,1:2),[],2)./sigma_pred,'k^','markerfacecol','k');
plot(smooth(max(sigmas(:,1:2),[],2)./sigma_pred,20,'rloess'),['-.k'],'linewid',2);


% plot(1:failureBegin-1,sigma_pred_ratio(1:failureBegin-1),'ks','markerfacecol','k','markersize',9);
% plot(failureBegin:n, sigma_pred_ratio(failureBegin:end),'ks','markersize',9);

if ~showPlatform
    scatter(1:n,sigma_pred_ratio,150,linspace(0,1,n),'fill','s');
else
    for pp = 1:length(platforms)
        ind = find(platform == platforms(pp));
        plot(ind,sigma_pred_ratio(ind),'s','color',platformColor{pp},'markerfacecol',platformColor{pp},'markersize',9);
    end
end


plot(smooth(sigma_pred_ratio,30,'rloess'),['-k'],'linewid',2);


ylim([0.3 3]); xlim([0 n+2]);
set(gca,'YScale','log','Ytick',[0.5 1 2:3]); 

xlabel('Session');
ylabel('Prediction ratio');
% ylabel('Prediction ratio = actual / predicted threshold');

Annotation(date,tf,glass);

SetFigure();


% %{
%% Correlations 

h = LinCorr(6,sigmas(:,1)./sigmas(:,2),sigma_pred_ratio,1,1,failureBegin);
xlabel('Log ( vesibular / visual threshold)');
ylabel('Log ( prediction ratio )');
SetFigure()

h = LinCorr(16,max(sigmas(:,1:2),[],2)./min(sigmas(:,1:2),[],2),sigma_pred_ratio,1,1,failureBegin);
xlabel('Mismatching = log ( max / min threshold)');
ylabel('Log ( prediction ratio )');
SetFigure()




%% 
% figure(7); clf

LinCorr(7,sigma_pred,sigmas(:,3),0,0,failureBegin);

% plot(sigma_pred(1:failureBegin-1),sigmas(1:failureBegin-1,3),'ko','markersize',10,'markerfacecol','k');
% plot(sigma_pred(failureBegin:end),sigmas(failureBegin:end,3),'ko','markersize',10,'linew',2);

xlabel('Predicted threshold)');
ylabel('Actual threshold');
plot([1 6],[1 6],'--k'); 

axis([1.1 6 1.1 6]); 
axis square;
SetFigure();

%%
LinCorr(17,sigma_pred,sigma_pred_ratio,1,1,failureBegin);

% plot(sigma_pred(1:failureBegin-1),sigmas(1:failureBegin-1,3),'ko','markersize',10,'markerfacecol','k');
% plot(sigma_pred(failureBegin:end),sigmas(failureBegin:end,3),'ko','markersize',10,'linew',2);

xlabel('Log (predicted threshold)');
ylabel('Log (prediction ratio)');
SetFigure();

%% 
LinCorr(8,sigmas(:,1),sigma_pred_ratio,1,1,failureBegin);

% figure(8); clf
% plot(sigmas(:,1),sigma_pred_ratio,'ok');
xlabel('Log (vestibular threshold)');
ylabel('Log (prediction ratio)');

SetFigure();

%%
LinCorr(9,sigmas(:,2),sigma_pred_ratio,1,1,failureBegin);

% figure(9); clf
% plot(sigmas(:,2),sigma_pred_ratio,'ok');
xlabel('Log (visual threshold)');
ylabel('Log (prediction ratio)');

SetFigure();

%%
LinCorr(10,min(sigmas(:,1:2),[],2),sigma_pred_ratio,1,1,failureBegin);

% plot(min(sigmas(:,1:2),[],2),sigma_pred_ratio,'ok');
xlabel('Log (Min threshold)');
ylabel('Log (prediction ratio)');

SetFigure();

%%
LinCorr(11,max(sigmas(:,1:2),[],2),sigma_pred_ratio,1,1,failureBegin);

% plot(max(sigmas(:,1:2),[],2),sigma_pred_ratio,'ok');
xlabel('Log (Max threshold)');
ylabel('Log (prediction ratio)');

SetFigure();

%%
keyboard;

%}


function Annotation(date,tf,glass)

% Same day
diff0 = find(diff(date)== 0);
diff0_begs = [diff0(1); diff0(1+find(diff(diff0)>1))];
diff0_ends = [diff0(diff(diff0)>1) ;diff0(end)]+1;
ylims = ylim;
xlims = xlim;

for j = 1:length(diff0_begs);
    plot([diff0_begs(j)-0.2 diff0_ends(j)+0.2],[ylims(1) ylims(1)],'k-','linewid',10);
end

% Not working
plot([70 70],ylims,'--k','linewi',2);
plot([74 74],ylims,'--k','linewi',2);


% Target first v.s. 3D glass
a1 = gca;
pos = get(gca,'Position');
a2 = axes('position',[pos(1) pos(2)+pos(4) pos(3) 0.03]); 

[xx, yy]=meshgrid(1:length(date),0); 
plot(xx(logical(tf)),yy(logical(tf)),'cs','markerfacecolor','c'); hold on;

[xx, yy]=meshgrid(1:length(date),1);
plot(xx(logical(glass)),yy(logical(glass)),'ms','markerfacecolor','m');

axis off;
xlim(xlims); ylim([0 1])
linkprop([a1 a2],'xlim');


function h = LinCorr(figN,x,y,logx,logy,seperate)
set(figure(figN),'position',[180 93 818 679]); clf;

nanInd = or(isnan(x),isnan(y));
x(nanInd) = [];
y(nanInd) = [];

if logx
    xx = log(x);
else 
    xx = x;
end

if logy
    yy = log(y);
else
    yy = y;
end

ind = [1 seperate length(x)];
markerfacecolors = {'k','none',[0.8 0.8 0.8]};
linestyles = {'b-','r-','k-'};

% All fitting
[r,p] = corr(xx,yy,'type','spearman');     
[para,S]=polyfit(xx,yy,1);

xxx = min(xx):0.1:max(xx);
Y = polyval(para,xxx);
h(10+length(ind)) = plot(xxx,Y,'k','linewidth',2);  hold on;

text(xxx(end),Y(end),sprintf('r = %3.3f \n p = %3.3f',r,p),'color','k');
scatter(xx,yy,70,linspace(0,1,length(xx)),'fill');


% Seperate draw points and fitting

for i = 1:(length(ind)-1)
    thisX = xx(ind(i):ind(i+1));
    thisY = yy(ind(i):ind(i+1));
    
%     h(i) =  scatter(thisX,thisY,100,linspace(0,1,length(thisX)),'fill');
    
%     h(i) = plot(thisX,thisY,'ko','markersize',7,'markerfacecol',markerfacecolors{i}); hold on;
    
    [r,p] = corr(thisX,thisY,'type','spearman');     
    [para,S]=polyfit(thisX,thisY,1);
    
    xxx = min(thisX):0.1:max(thisX);
    Y = polyval(para,xxx);
    h(i+10) = plot(xxx,Y,[linestyles{i}],'linewidth',2);

    text(xxx(end),Y(end),sprintf('r = %3.3f \n p = %3.3f',r,p),'color',linestyles{i}(1));
end







