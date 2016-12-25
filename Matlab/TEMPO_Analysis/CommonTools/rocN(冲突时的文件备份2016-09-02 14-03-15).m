function [a,bestz] = rocN(x,y,N)
% ROC	computes area under ROC given distributions x and y
%	uses N points to construct the ROC
%	usage area = rocN(x,y,N)
%
%   x, +; y, -
%       Warning -- it gives a roc of NaN if there are any Nans!
% 
if nargin < 3
  N = 100;
end

% [m n] = size(x);
% x = reshape(x,1,m*n);
% [m n] = size(y);
% y = reshape(y,1,m*n);
x = x(:);
y = y(:);

if sum(isnan(x)) + sum(isnan(y)) > 0 || isempty(x) || isempty(y) % Any NaN or any []
    a = NaN;
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

a = trapz(fa,hit);
%}


% %{  
% HH20141001
z = sort([-inf ; x; y ; inf]); % Meaningful threshold steps

if length(z)<N   % No need to step that much steps
    N = length(z);
else
    zlo = min([min(x(:)) min(y(:))]);
    zhi = max([max(x(:)) max(y(:))]);

    z = [-inf linspace(zlo,zhi,N) inf];
end

fa = zeros(1,N); hit = fa; % Preallocation

for i = 1:N
   fa(N-i+1) = sum(y > z(i));
   hit(N-i+1) = sum(x > z(i));
end

fa = fa/length(y);
hit = hit/length(x);

a = trapz(fa,hit);

if nargout == 2 % Return the best threshold. @HH20150204
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
subplot(1,2,1);
plot(fa,hit,'.-r'),axis('square'),xlabel('FA'),ylabel('Hit'); hold on;
plot([0 1],[0 1],'k--');

subplot(1,2,2);

xx = linspace(min([x;y]),max([x;y]),20);
[posCont,posX] = hist(x,xx);
[negCont,negX] = hist(y,xx);

bar(xx,posCont,'g'); hold on; 
bar(xx,-negCont,'r');

maxabsY = max(abs(ylim));
ylim(maxabsY * [-1 1]);

if exist('bestz') & ~isnan(bestz)
    plot([bestz bestz],ylim,'k--','linew',2);
    subplot(1,2,1); plot(fa(ind),hit(ind),'ob');
end
% pause

%}