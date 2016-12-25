function [ret]=pois(X,a)
A=a(1);
mu=a(2);
ret=A*(mu.^X.*exp(-mu))./gamma(X+1);
