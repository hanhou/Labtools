%% original code
clear;

%      duration   amplitude   num_sigs   notes
paras = {
%     2   0.13  6    'Gu 2006 2010 tuning (combined)'
%           2   0.2   4.5  'Gu 2012 heading (combined)'
%           1   0.2   2    'HeTao heading (visual only)'
%           2 0.18 6 'Polo training'
%           1.5   0.2  3.5   'Polo heading (combined)'
          1.5   0.1  4  'Polo heading (combined)'
          1.5   0.2  3.5  'Polo heading (combined)'
          1.5   0.14  3.8  'Minimoog'
          };


figure(100); clf;
set(100,'position',[100 100 900 700]); 

maxDuration =  max(cell2mat(paras(:,1)));

for i=1:size(paras,1)
    
    
    duration = paras{i,1};
    step=0.005;
    t = 0:step:duration;
    
    ampl = paras{i,2};
    num_sigs = paras{i,3};
    
    %pos = ampl*0.5*(erf(2*num_sigs1/3*(t-1)) + 1);
    pos=ampl*0.5*(erf(sqrt(2)*num_sigs*(t-duration/2)/duration) + 1);
    %pos = ampl*normcdf(t,1,1/num_sigs);
    
    subplot(size(paras,1),4,1+4*(i-1));
    plot(t,pos,'linewidth',2);
    ylim([min(pos) max(pos)]);
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
    title(sprintf('Max pos = %6.2f m', max(pos))); 
    
    veloc = diff(pos)/step;
    subplot(size(paras,1),4,2+4*(i-1));
    plot(t(1:length(t)-1),veloc,'r-','linewidth',2); 
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
    ylim([min(veloc) max(veloc)]);
    title(sprintf('Max veloc = %6.2f m/s', max(veloc))); 
    max(veloc)
 
    
    accel = diff(veloc)/step/9.8;
    
    subplot(size(paras,1),4,3+4*(i-1));
    plot(t(1:length(t)-2),accel,'g-','linewidth',2); 
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
%     ylim([min(accel) max(accel)]);
    ylim([-0.12 0.12]);
    title(sprintf('Max accel = %6.2f g', max(accel)));
    max(accel)
    
    
    subplot(size(paras,1),4,4+4*(i-1));
    axis('off')
    buff = sprintf('%s\n\nduration = %6.2f\nampl = %6.2f\nnum sigs = %6.2f', paras{i,4}, duration, ampl, num_sigs);
    text(.1,.8,buff);
end

%%
% %{
am = 0.04:0.01:0.32;
sig = 3:0.1:8;
% bb(:,1) = am';

duration = 1.5;
step=0.005;
t = 0:step:duration;


for i = 1:length(am)
    for j = 1:length(sig)

        ampl = am(i);
        num_sigs = sig(j);

        % pos = ampl*0.5*(erf(2*num_sigs/3*(t-1)) + 1);
        pos=ampl*0.5*(erf(sqrt(2)*num_sigs*(t-duration/2)/duration) + 1);

        veloc = diff(pos)/step;
        accel = diff(veloc)/step/9.8;

        aa(i,j) = max(accel);
%         bb(i,j+1) = max(accel);
    end
end
figure(3); 
contourf(sig,am,aa); 
title(sprintf('duration = %3.2f, accel',duration));

% bb
colorbar;
xlabel('sigma');
ylabel('amol(m)');

%}

