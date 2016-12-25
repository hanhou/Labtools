aa=dlmread('eyetracechoice.txt');
h{1}='r-';h{2}='b-';h{3}='g-';h{4}='y-';h{5}='k-';h{6}='y--';h{7}='g--';h{8}='b--';h{9}='r--';
dim = size(aa);
%for i=1:9
%     pp(i,:) = median(aa(i:9:end,1:400));
%     plot( median(aa(i:9:end,1:400)),h{i},'LineWidth',2 );
%     xlim([0 400]);
%     ylim([-0.1,0.1]);
%     hold on;
%end

% plot(pp','LineWidth',2);
stim = 1;
% subplot(2,2,1);
% plot( (aa(stim:3:end,1:400)-aa(stim:3:end,401:800))' );
% hold on;
% plot( median(aa(stim:3:end,1:400)-aa(stim:3:end,401:800)),'r-','LineWidth',3 );
% xlim([0 400]);
% ylim([-0.2,0.2]);
% xlabel('Left Choice');
% subplot(2,2,2);
% plot( (aa(stim:3:end,801:1200)-aa(stim:3:end,1201:1600))');
% hold on;
% plot( median(aa(stim:3:end,801:1200)-aa(stim:3:end,1201:1600)),'b-','LineWidth',3 );
% hold on;
% xlim([0 400]);
% ylim([-0.2,0.2]);
% xlabel('Right Choice');
% subplot(2,2,3);
% plot( median(aa(stim:3:end,1:400)-aa(stim:3:end,401:800)),'r-','LineWidth',1 );
% hold on;
% plot( median(aa(stim:3:end,801:1200)-aa(stim:3:end,1201:1600)),'b-','LineWidth',1 );
% hold on;
% xlim([0 400]);
% ylim([-0.2,0.2]);
xxlim=0.1;
yy = 1:5:2000
LL = median(aa(1:3:end,1:400));
RL = median(aa(1:3:end,401:800));
LR = median(aa(1:3:end,801:1200));
RR = median(aa(1:3:end,1201:1600));
subplot(1,3,1);
plot( LL,yy, 'r-', 'LineWidth', 2 );
hold on;
plot( RL,yy, 'r-', 'LineWidth', 1 );
hold on;
plot( LR,yy, 'b-' , 'LineWidth', 2);
hold on;
plot( RR,yy, 'b-', 'LineWidth', 1);
hold on;
xlim([-xxlim,xxlim]);
xlabel('Horizontal eye position');
ylabel('stimulus duration');
title('ves');

LL = median(aa(2:3:end,1:400));
RL = median(aa(2:3:end,401:800));
LR = median(aa(2:3:end,801:1200));
RR = median(aa(2:3:end,1201:1600));
subplot(1,3,2);
plot( LL,yy, 'r-', 'LineWidth', 2 );
hold on;
plot( RL,yy, 'r-', 'LineWidth', 1 );
hold on;
plot( LR,yy, 'b-' , 'LineWidth', 2);
hold on;
plot( RR,yy, 'b-', 'LineWidth', 1);
hold on;
xlim([-xxlim,xxlim]);
xlabel('Horizontal eye position');
ylabel('stimulus duration');
title('vis');

LL = median(aa(3:3:end,1:400));
RL = median(aa(3:3:end,401:800));
LR = median(aa(3:3:end,801:1200));
RR = median(aa(3:3:end,1201:1600));
subplot(1,3,3);
plot( LL,yy, 'r-', 'LineWidth', 2 );
hold on;
plot( RL,yy, 'r-', 'LineWidth', 1 );
hold on;
plot( LR,yy, 'b-' , 'LineWidth', 2);
hold on;
plot( RR,yy, 'b-', 'LineWidth', 1);
hold on;
xlim([-xxlim,xxlim]);
xlabel('Horizontal eye position');
ylabel('stimulus duration');
title('comb');





