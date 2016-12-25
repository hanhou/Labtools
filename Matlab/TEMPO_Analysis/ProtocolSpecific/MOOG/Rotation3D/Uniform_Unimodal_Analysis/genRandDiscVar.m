function X=genRandDiscVar(freq, n, s)
% function X=genRandDiscVar(freq, n, s)
% generate discrete random variables
%
% freq: a vector of length m that contains frequencies of m discrete variate values
% n: the number of samples to generate
% s: state for random number generator reset

DEBUG=0;

if s>0
    rand('state',s);
end

m=length(freq);
cum=cumsum(freq);
cumlativeFreq=cum/cum(end);

X=rand([1 n]);
for i=1:n
    while X(i)==1
        X(i)=rand;
    end
end

for i=1:n
    for j=1:m
        if X(i)<=cumlativeFreq(j);
            X(i)=j;
            break;
        end
    end
end

if DEBUG
    figure(1);
    xh=1:m;
    yh=hist(X,xh);

    v=zeros(2,m);
    v(1,:)=yh/max(yh);
    v(2,:)=freq/max(freq);

    bar(v','group');
    axis([0 length(xh)+1 0 1]);
end