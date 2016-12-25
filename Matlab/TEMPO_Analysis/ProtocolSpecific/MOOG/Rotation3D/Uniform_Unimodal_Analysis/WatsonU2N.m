function U2N=WatsonU2N(v)
% function U2N=WatsonU2N(v)
% computes Watson's U2 statistic for goodness of fit tests on a circle (Watson, 1961)
%
% v: a vector of length N, and v(i)=F(x(i)) where F(x) is a continuous
%    cummulative distribution function for a population (theoretical)
%    evaluated at sample points x.
%
% return
% U2N: Watson's U2 goodness-of-fit statistic
%      if the sample U2 exceeds the critical value, H0 is rejected, i.e. not a good fit
%      otherwise, H0 cannot be rejected, i.e. a good fit
%
% tested using Example 4.10.1 in Batschelet (1981)

N=length(v);
v=sort(v);
vm=mean(v);
i=1:N;
c=2*i-1;
v1=c*v'/N;
v2=v*v';
U2N=v2-v1+N/3-N*(vm-0.5)*(vm-0.5);
