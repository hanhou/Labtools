function f=compKDEvonMises(xf,xh,k)
%compute a kernel density estimate
m=length(xh);
f=zeros(size(xf));
for i=1:m
    v=vonMises(xf,[xh(i) k]);
    f=f+v;
end
%h=1/sqrt(k); %this is Fisher's method
h=acos(1-1/(2*k)); %this is my method
f=f/h/m;
