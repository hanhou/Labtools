load z:\Users\Jacob\trajectory.mat

i = 14; j=10;
a = [zeros(1,i) predLateral -.001:.00011:0];
a = interp(a,500);
a = decimate(a,127+i+j);
a = a.*100;
b = [diff(a) 0];
b = b.*200;
b = boxcarfilter(b,40);
c = [diff(b) 0];
c = c.*100
figure(3);
clf;
hold on
plot(a,'r')
figure(4);
clf;
hold on
plot(b,'r')
figure(5);
clf;
hold on
plot(c,'r')