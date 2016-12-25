figure(1);
x=0:1:360;
 y=-90:1:90;yy=y/180*pi;
contourf(x,sin(yy),example_fit);
set(gca,'ydir','reverse');
caxis([30, 170]);
colorbar;
axis off;hold on;
xp=[109];yp=[-3];yyp=yp/180*pi;
plot(xp,yyp,'w+');

figure(2);
xx=1:1:360;
contourf(xx,sin(yy),example_fisher);
set(gca,'ydir','reverse');
colorbar;
axis off;