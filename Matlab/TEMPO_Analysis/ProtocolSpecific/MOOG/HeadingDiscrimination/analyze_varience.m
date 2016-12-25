% % % for 100% visual
% aa=dlmread('vis.txt');
% dim=size(aa);
% 
% % re-arragne data to avoid NaN value
% for i = 1 : 61
%     for j = 1 : 11
%         rr(j+11*(i-1), 1) = aa(i,j); % for sets with 11 data points
%         rr(j+11*(i-1), 2) = aa(i,j+11);
%     end
% end
% for i = 1 : 41
%     for j = 1 : 9
%         rr(j+9*(i-1)+61, 1) = aa(i+61,j); % for sets with 11 data points
%         rr(j+9*(i-1)+61, 2) = aa(i+61,j+11);
%     end
% end
% 
% figure;
% %plot( aa(:,1:dim(2)/2), aa(:,dim(2)/2+1 : end), 'r.');
% plot( rr(:,1), rr(:,2), 'r.');
% xlim([0,150]);
% ylim([0,150]);
% title('100% Coherence Visual');
% hold on;
% plot( [0,150], [0,150], 'k--');
% 
% pp = polyfit(rr(:,1), rr(:,2), 1);
% xx = 0:1:150;
% yy = pp(1)*xx + pp(2);
% plot(xx,yy,'r-');
% text(120,70,num2str(pp(1)) );
% xlabel('mean');
% ylabel('variance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for 3 conditions
aa=dlmread('all.txt');
ss=dlmread('sign.txt');
dim=size(aa);

figure;
plot( log10(aa(:,1:9)), log10(aa(:,10 : 18)), 'b.');
hold on;
plot( log10(aa(:,19:27)), log10(aa(:,28 : 36)), 'r.');
plot( log10(aa(:,37:45)), log10(aa(:,46 : 54)), 'g.');
xlim([-1,3]);
ylim([-1,3]);
title('3 conditions');
plot( [-1,3], [-1,3], 'k-');

pp1 = polyfit(log10(aa(:,1:9)), log10(aa(:,10 : 18)), 1);
pp2 = polyfit(log10(aa(:,19:27)), log10(aa(:,28 : 36)), 1);
pp3 = polyfit(log10(aa(:,37:45)), log10(aa(:,46 : 54)), 1);
xx = -1:0.1:3;
yy1 = pp1(1)*xx + pp1(2);
yy2 = pp2(1)*xx + pp2(2);
yy3 = pp3(1)*xx + pp3(2);
plot(xx,yy1,'b--','LineWidth',2);
plot(xx,yy2,'r--','LineWidth',2);
plot(xx,yy3,'g--','LineWidth',2);
text(-0.5,2.5,num2str(pp1(1)) );
text(-0.5,2.3,num2str(pp2(1)) );
text(-0.5,2.1,num2str(pp3(1)) );
xlabel('mean');
ylabel('variance');

% %%----------------------------------------------------------------
% % linear
% figure;
% plot( aa(:,1:9), aa(:,10 : 18), 'b.');
% hold on;
% plot( aa(:,19:27), aa(:,28 : 36), 'r.');
% plot( aa(:,37:45), aa(:,46 : 54), 'g.');
% xlim([0,150]);
% ylim([0,150]);
% title('3 conditions');
% plot( [0,150], [0,150], 'k--');
% 
% pp1 = polyfit(aa(:,1:9), aa(:,10 : 18), 1);
% pp2 = polyfit(aa(:,19:27), aa(:,28 : 36), 1);
% pp3 = polyfit(aa(:,37:45), aa(:,46 : 54), 1);
% xx = 0:1:150;
% yy1 = pp1(1)*xx + pp1(2);
% yy2 = pp2(1)*xx + pp2(2);
% yy3 = pp3(1)*xx + pp3(2);
% plot(xx,yy1,'b-','LineWidth',2);
% plot(xx,yy2,'r-','LineWidth',2);
% plot(xx,yy3,'g-','LineWidth',2);
% text(100,50,num2str(pp1(1)) );
% text(100,45,num2str(pp2(1)) );
% text(100,40,num2str(pp3(1)) );
% xlabel('mean');
% ylabel('variance');
