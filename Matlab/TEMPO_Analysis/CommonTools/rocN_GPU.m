function a = rocN_GPU(x,y,N)
% ROC	computes area under ROC given distributions x and y
%	uses N points to construct the ROC
%	usage area = rocN(x,y,N)
%
%   x, +; y, -
%       Warning -- it gives a roc of NaN if there are any Nans!

if nargin < 3
  N = gpuArray(100);
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
    N = gpuArray(length(z));
else
    zlo = min([min(x(:)) min(y(:))]);
    zhi = max([max(x(:)) max(y(:))]);

    z = [-inf gpuArray.linspace(zlo,zhi,N) inf];
end

fa = gpuArray.zeros(1,N); 
hit = fa; % Preallocation

i = 1;
while i <= N 
   fa(N-i+1) = sum(y > z(i));
   hit(N-i+1) = sum(x > z(i));
   i = i + 1;
end

fa = fa/length(y);
hit = hit/length(x);

a = trapz(fa,hit);  

%}

% uncomment next 2 lines if you want to see the plot
% figure(1); hold on;
% plot(fa,hit,'.-r'),axis('square'),xlabel('FA'),ylabel('Hit');

