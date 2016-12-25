a = [ 17.5 16.6 11.5 11.5 8.9 8.4 6.75 5.4 4.5 ];
b = [ 10.85 9.85 7.95 6.2 6.05 5.9 9.4 9.5 13.3 ];
a = [36.5 33.35 28.35 24.25 17.8 14.65 12.95 13.1 18.5];

numparamq = 2;
numparaml = 1;

range = [-2:.5:2];

pf = polyfit(range,a,numparamq);
q = polyval(pf,range);
pf1 = polyfit(range,a,numparaml);
l = polyval(pf1,range);

sseq = sum((a - q).^2);
ssel = sum((a - l).^2);

F = ((ssel-sseq)/(numparamq-numparaml))/(sseq/(9-numparamq))

pvalue = 1 - fcdf(F,(numparamq-numparaml),(9-numparamq))

figure(2);
clf;
hold on;
plot(a);
plot(q,'r');
plot(l,'g');