aa = dlmread('shiftwindow.txt'); 
% 1 :11  threshold 
% 12:22  CP
% 23:33  p
dim = size(aa);
maxpsth = dlmread('maxpsth.txt'); 
minpsth = dlmread('minpsth.txt');
% 21 : 60 two seconds stimulus duration

% normalize threshold first
for i = 1 : dim(1)
    % normalize threshold
    aa_norm(i,1:11) = aa(i,1:11) / max(aa(i,1:11));
    [rr,pp] = corrcoef(aa_norm(i,1:11), aa(i,12:22) );
    r(i) = rr(1,2);
    p(i) = pp(1,2);
    % sort neuron by where is the smallest threshold
    mm = find( aa_norm(i,1:11) == min(aa_norm(i,1:11)) ); 
    thresh_min_index(i) = mm(1);
    CP_corespond_minthreshold(i) = aa(i,mm(1)+11);
    % sort neuron by where is the largest CP
    cc = find( aa(i,12:22) == max(aa(i,12:22)) );
    CP_max_index(i) = cc(1);  % in case two same minium values
    ccmin = find( aa(i,12:22) == min(aa(i,12:22)) );  % also find smallest CP
    CP_min_index(i) = ccmin(1);  % in case two same minium values
end
aa(:,1:11) = aa_norm(:,1:11);

% sort neuron according to smallest threshold
[ss, ii] = sort(thresh_min_index);
for i = 1 : dim(1)
    aa_sort(i,:) = aa(ii(i),:);
    maxpsth_sort(i,:) = maxpsth(ii(i),:);
    minpsth_sort(i,:) = minpsth(ii(i),:);
end

% % % plot contour map, row is cell number, column is threshold or CP
% figure(2);
% subplot(1,2,1);
% contourf(aa_sort(:,1:11));
% colorbar;
% xlabel('shiftwindow');
% set(gca,'xtick', [1,2,3,4,5,6,7,8,9,10,11]);
% ylabel('cell #');
% title('Threshold');
% 
% subplot(1,2,2);
% contourf(aa_sort(:,12:22));
% colorbar;
% set(gca,'xtick', [1,2,3,4,5,6,7,8,9,10,11]);
% xlabel('shiftwindow');
% ylabel('cell #');
% title('CP');

% % % plot contour map, row is cell number, column is threshold or psth
% figure(2);
% subplot(1,2,1);
% contourf(aa_sort(:,1:11));
% colorbar;
% xlabel('shiftwindow');
% set(gca,'xtick', [1,2,3,4,5,6,7,8,9,10,11]);
% ylabel('cell #');
% title('Threshold');
% 
% subplot(1,2,2);
% contourf(minpsth_sort(:,21:60));
% colorbar;
% %set(gca,'xtick', []);
% xlabel('stimulus duration');
% ylabel('cell #');
% title('PSTH mindirection');

% subplot(1,2,2);
% contourf(maxpsth_sort(:,21:60));
% colorbar;
% %set(gca,'xtick', []);
% xlabel('stimulus duration');
% ylabel('cell #');
% title('PSTH maxdirection');
