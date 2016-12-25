a=rand(10,10)
[n,bin]=histc(a,0:0.5:1); 

col = 3; group = 2;
x = a(bin(:,col) == group , col)