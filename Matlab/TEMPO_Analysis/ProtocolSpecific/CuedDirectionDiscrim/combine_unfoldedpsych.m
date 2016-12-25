datadir = 'Z:\\Data\Tempo\\Baskin\Analysis\Unfolded_CuedDirecDiscrim_Psych\';
d = ls(datadir);
d = d(3:end-1,:);
coher = [-16 -8 -4 -2 0 0 2 4 8 16];
cum_psychdata = [];
for i = 1:length(d)
    c = importdata(sprintf('%s%s',datadir,d(i,:)));
    cum_psychdata(:,:,end+1) = c(2:4,1:10); %indices of pct correct at each validity
end
plot_colors = {'bs','ro','g>'};
plot_lines = {'b','r','g'};
mean_psychdata = mean(cum_psychdata,3);
std_psychdata = std(cum_psychdata,0,3);
figure; hold on;
for i = 1:3
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_psychdata(i,:);
    fit_data(:,3) = size(cum_psychdata,3);
    [group_alpha(i) group_beta(i)] = logistic_fit(fit_data);
    fit_x = [min(coher):0.1:max(coher)];
    fit_y(i,:) = logistic_curve(fit_x,[group_alpha(i) group_beta(i)]);
    errorbar(coher,mean_psychdata(i,:),2.*std_psychdata(i,:)./sqrt(length(d)),plot_colors{i});
    plot(fit_x,fit_y(i,:),plot_lines{i});
end
% temp = [coher' mean_psychdata'; coher' std_psychdata'; fit_x' fit_y'];
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_psych.txt', '-ascii', 'temp')

yneu = mean_psychdata(i,:);
neuvar = yneu.*(1-yneu);
v = 2/3; vvar = v*(1-v);
valcue = [0 0 0 0 0 1 1 1 1 1]; 
invcue = 1-valcue;
both = [1-v 1-v 1-v 1-v 1-v v v v v v];
yinv = invcue.*(neuvar./(vvar+neuvar)) + yneu.*(vvar./(vvar+neuvar));
yval = valcue.*(neuvar./(vvar+neuvar)) + yneu.*(vvar./(vvar+neuvar));
yboth = both.*(neuvar./(vvar+neuvar)) + yneu.*(vvar./(vvar+neuvar));
plot(coher,yinv,'bs', 'MarkerFaceColor',[0 0 1]);
plot(coher,yval,'g>', 'MarkerFaceColor',[0 1 0]);
axis('tight');