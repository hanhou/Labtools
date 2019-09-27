%% Calculate Gaussian profiles based on the method used in MoogDots.
% HH 2012
 
clear;

%      duration   amplitude   num_sigs   notes
paras = {
          2   0.13  6    'Gu 2006 2010 tuning (combined)'
%           2   0.2   4.5  'Gu 2012 heading (combined)'
%           1   0.2   2    'HeTao heading (visual only)'
%           2 0.18 6 'Polo training'
          1.5   0.2  3.5   'Polo heading (original)'
%           1.5   0.1  4  'Polo heading (combined)'
%           2.0   0.2  4.0  'Polo heading (combined)'
%           1.5   0.14  3.8  'Minimoog'
          1.5   0.1  4.9  'Change sigma 1' 
          1.5   0.25  2.3  'Change sigma 2' 
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
    
    %pos = ampl*0.5*(erf(2*num_sigs1/3*(t-1)) + 1); % Erroneous
    pos=ampl*0.5*(erf((t-duration/2)/((duration/2/num_sigs)*sqrt(2))) + 1);  % Note that sigma = duration/2/num_of_sigma !!
%     pos = ampl * normcdf(t,duration/2,duration/2/num_sigs); % This is much simplier
    
    subplot(size(paras,1),4,1+4*(i-1));
    plot(t,pos,'linewidth',2);
    ylim([min(pos) max(pos)]);
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
    title(sprintf('Max pos = %6.2f m', max(pos))); 
    
    veloc = diff(pos)/step;
    subplot(size(paras,1),4,2+4*(i-1));
    plot(t(1:length(t)-1),veloc,'r-','linewidth',2); 
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
    ylim([0 max(veloc)]);
    title(sprintf('Max veloc = %6.2f m/s', max(veloc))); 
    t_max_veloc = t(veloc == max(veloc)); 
 
    
    accel = diff(veloc)/step/9.8;
    
    subplot(size(paras,1),4,3+4*(i-1));
    plot(t(1:length(t)-2),accel,'g-','linewidth',2); 
    xlim([duration/2-maxDuration/2 duration/2+maxDuration/2]);
%     ylim([min(accel) max(accel)]);
    ylim([-0.12 0.12]);
    title(sprintf('Max accel = %6.2f g', max(accel)));
    t_max_accel = t(accel == max(accel));
    
    
    subplot(size(paras,1),4,4+4*(i-1));
    axis('off')
    buff = sprintf('%s\n\nduration = %6.2f\nampl = %6.2f\nnum sigs = %6.2f\nsig = %g ms\ndt(v-a peak) = %g ms', ...
                    paras{i,4}, duration, ampl, num_sigs, duration/2/num_sigs, (t_max_veloc(1)-t_max_accel(1))*1000);
    text(.1,.8,buff);
end


%% Parameter selection for sigma change

duration = 1.5;
ampls = 0.1:0.03:0.3;
num_sigss = 1.0:0.1:7.0;

step=0.002;
t = 0:step:duration;

max_accels = []; dt_avpeaks = [];

figure(200);
col = [0.9 0.9 0.9];

for i = 1:length(ampls)
    for j = 1:length(num_sigss)
        
        ampl = ampls(i);
        num_sigs = num_sigss(j);
        
        %pos = ampl*0.5*(erf(2*num_sigs1/3*(t-1)) + 1); % Erroneous
        pos=ampl*0.5*(erf((t-duration/2)/((duration/2/num_sigs)*sqrt(2))) + 1);  % Note that sigma = duration/2/num_of_sigma !!
        %     pos = ampl * normcdf(t,duration/2,duration/2/num_sigs); % This is much simplier
        
        veloc = diff(pos)/step;
        t_max_veloc = t(veloc == max(veloc));
        
        accel = diff(veloc)/step/9.8;
        t_max_accel = t(accel == max(accel));
        
        max_accels(i,j) = max(accel);
        dt_avpeaks(i,j) = 1000 * (t_max_veloc(1) - t_max_accel(1));
        
    end

    col = col - [0.07 0.07 0];
    plot(max_accels(i,:), dt_avpeaks(i,:),'color',col,'linew',2); hold on;
    
end
title(['dur = ' num2str(duration)]);
legend(num2str(ampls'));
xlabel('Max acceleration (g)'); ylabel('dt va peak'); SetFigure(15);
xlim([0 0.15]); ylim([100 500]);

%% Original code from Gu (which is erroneous)
%{
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

%% ZQH (this makes the amp accurate even when num_of_sigma is very small)
%{
% MOOG parameters
figure(111);clf
duration = 1.5; num_of_sigma = 1; Amp = 0.14;

% function parameters
miu = duration/2; sigma = duration/2/num_of_sigma; 
syms t;
A = Amp/int(exp(-(t-miu).^2/(2*sigma^2)),0,duration)

% functions
syms t;
v=A*exp(-(t-miu).^2/(2*sigma^2));
a=diff(v);
s=int(v,0,t);

% plot curves
t=0:0.01:duration;
plot(t,subs(s),'r','linewidth',3);hold on;
plot(t,subs(v),'b','linewidth',3);
plot(t,subs(a),'g','linewidth',3);
%}
