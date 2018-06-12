function [auROC,bestz,perm] = rocN(x,y,N,permuteN)
% ROC	computes area under ROC given distributions x and y
%	uses N points to construct the ROC
%	usage area = rocN(x,y,N)
%
%   x, +; y, -
%       Warning -- it gives a roc of NaN if there are any Nans!
% 
if nargin < 3 || isempty(N)
    N = 100;
end

if nargin < 4
    permuteN = 0;
end


% [m n] = size(x);
% x = reshape(x,1,m*n);
% [m n] = size(y);
% y = reshape(y,1,m*n);
x = x(:);
y = y(:);

if sum(isnan(x)) + sum(isnan(y)) > 0 || isempty(x) || isempty(y) % Any NaN or any []
    auROC = NaN;
    bestz = NaN;
    perm.pValue = NaN;
    perm.std = NaN;
    perm.auROC = NaN;
    return;
end


%{
zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);

z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);

for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end

fa = fa/length(y);
hit = hit/length(x);

fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;

a1 = trapz(fa,hit);
%}


% %{  
% HH20141001
% Modified by HH20160623

z = sort([-inf ; x; y ; inf]); % Meaningful threshold steps

if length(z)<N   % No need to step that many steps
    N = length(z);
else  
%     zlo = min([min(x(:)) min(y(:))]);
%     zhi = max([max(x(:)) max(y(:))]);
    z = [-inf linspace(z(2),z(end-1),N-2) inf]; % Redistribute Zs between min(x,y) and max(x,y) while keeping the length of z to be N. HH20160623
end

fa = zeros(1,N); hit = fa; % Preallocation

for i = 1:N  
   fa(N-i+1) = sum(y > z(i));
   hit(N-i+1) = sum(x > z(i));
end

fa = fa/length(y);
hit = hit/length(x);

auROC = trapz(fa,hit);

% Permutation here using matrix form. HH20180607
if permuteN > 0
    
    % Generate randPermX and Y
    [~,randPermSortIndex] = sort(rand(length([x;y]),permuteN),1); % This is really brilliant.
    tempXY = [x;y];
    randPermXY = tempXY(randPermSortIndex);
    randPermX = randPermXY(1:length(x),:);
    randPermY = randPermXY(length(x)+1:end,:);
    
    % Do ROC in batch
    faPerm = zeros(N,permuteN); hitPerm = faPerm; % Preallocation
    
    for i = 1:N
        faPerm(N-i+1,:) = sum(randPermY > z(i),1);
        hitPerm(N-i+1,:) = sum(randPermX > z(i),1);
    end
    
    faPerm = faPerm/length(y);
    hitPerm = hitPerm/length(x);
    
    perm.auROCPerm = nan(permuteN,1);
    for n = 1:permuteN
        perm.auROCPerm(n) = trapz(faPerm(:,n),hitPerm(:,n));
    end
    
    perm.pValue = sum(perm.auROCPerm > auROC) / permuteN; % One tail
    if auROC < 0.5, perm.pValue = 1 - perm.pValue; end
    perm.pValue = perm.pValue * 2; % Two tail
    
    perm.std = std(perm.auROCPerm);
else
    perm.pValue = NaN;
    perm.std = NaN;
    perm.auROCPerm = NaN;
end

if nargout > 1 % Return the best threshold. @HH20150204
    [correct_rate,ind] = max(hit - fa);
    if correct_rate > 0.1 
        bestz = z(ind);
    else
        bestz = nan;
    end
end
%}

%% Uncomment the following lines if you want to see the plot
%{
figure(1);  clf
subplot(1,3,1);
plot(faPerm,hitPerm,'color',0.8*ones(1,3)); hold on;
plot(fa,hit,'.-r'),axis('square'),xlabel('FA'),ylabel('Hit');
plot([0 1],[0 1],'k--');

subplot(1,3,2);
hist(perm.auROCPerm,20); hold on
plot([auROC auROC],ylim,'r','linew',2)
xlim([0 1])
title(sprintf('p = %g',perm.pValue));

subplot(1,3,3);

xx = linspace(min([x;y]),max([x;y]),20);
[posCont,posX] = hist(x,xx);
[negCont,negX] = hist(y,xx);

bar(xx,posCont,'g'); hold on; 
bar(xx,-negCont,'r');

maxabsY = max(abs(ylim));
ylim(maxabsY * [-1 1]);

if exist('bestz') & ~isnan(bestz)
    plot([bestz bestz],ylim,'k--','linew',2);
    subplot(1,3,1); plot(fa(ind),hit(ind),'ro');
end
% pause

%}