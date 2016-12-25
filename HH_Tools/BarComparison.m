function result = BarComparison(data,varargin)
% Draw bar comparison, recursive (paired) t-test, etc.
% Input: data   column for varibles, row for observations
%
% 20160607 HH

% ------ Parse input parameters -------
paras = inputParser;

addOptional(paras,'figN',555);
addOptional(paras,'axes',[]);

addOptional(paras,'Colors',{'b','r','g'}); % (originally for comparing PSTH of LIP neurons)
addOptional(paras,'LineStyles',{'-'});

parse(paras,varargin{:});

Colors = paras.Results.Colors;

% --------- End input parser ----------


% Basic statistics
result.means = mean(data,1);
result.stds = std(data,[],1);
result.sems = result.stds/sqrt(size(data,1));
n_var = size(data,2);

% Draw bars
figure(paras.Results.figN); clf; hold on; xlim([0.5 n_var+0.5]);

for i = 1:n_var
    col_ind = 1+mod(i-1,length(Colors));
    
    bar(i,result.means(i),0.7,'facecol',Colors{col_ind},'edgecol','none'); 
    herr = errorbar(i,result.means(i),result.sems(i),'linestyle','none','linewidth',3,'col',Colors{col_ind});
    errorbar_tick(herr,13);
end
        
% Overdraw raw data
plot(1:n_var,data,'ko-');

% Recursive paired t-test
ps_ttest = nan(size(data,2));
ps_vartest = nan(size(data,2));
name = 1:size(data,2);

for g1 = 1:size(data,2)
    for g2 = g1 + 1:size(data,2)
        % Paired ttest
        [~,p] = ttest(data(:,g1),data(:,g2));
        ps_ttest(g1,g2) = p;
        ps_ttest(g2,g1) = p;
        
        % Var test
        [~,p_var] = vartest2(data(:,g1),data(:,g2));
        ps_vartest(g1,g2) = p_var;
        ps_vartest(g2,g1) = p_var;
    end
end
%  sum(sum(ps<0.05))

result.ps_ttest = [[nan;name'] [name;ps_ttest] ];
result.ps_vartest = [[nan;name'] [name;ps_vartest] ];



