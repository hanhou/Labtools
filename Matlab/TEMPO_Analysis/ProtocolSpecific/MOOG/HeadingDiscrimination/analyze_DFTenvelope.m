aa = dlmread('all.txt');
dim = size(aa);

t = 0:0.05:2;
ampl = 0.30;
num_sigs = 4;
pos = ampl*0.5*(erf(2*num_sigs/3*(t-1)) + 1);
veloc = diff(pos)/0.01;

figure(2);
set(2,'Position', [5,5 1000,680], 'Name', 'Envelope');
orient landscape;
axis off;
x_time = 1:46;
x_start = [4,4];
x_stop = [43,43];
y_marker = [0,1];

% for k = 1:182
%     HT = aa(k,104:end-3) - mean(aa(k,104:105));% remove the mean to analyze just the peak or else the center of mass will be biased
%     HT(HT<0) = 0;  
%     H_time = 1:.05:3;
%     for hh = 1:length(HT)
%         CMsum(hh) = HT(hh)*H_time(hh);
%     end
%     if sum(HT) > 0
%         center_mass = sum(CMsum)/sum(HT);
%         peak_t(k) = find( abs(H_time - center_mass) == min(abs(H_time - center_mass)) );% point which is closest to the center of mass is used.
%     end
% end

for j = 1 : 4     
    for i = 1 : 4      
        axes('position',[0.24*(j-1)+0.05 (0.98-0.22*i) 0.2 0.18]);  
        % raw data
        max_norm = max( aa((j-1)*4+i, 22:67) );
        bar( x_time, aa((j-1)*4+i, 22:67)/max_norm,'w' );  
        hold on;     
        % hilbert
        plot(x_time, (aa((j-1)*4+i, 101:end)+mean(aa((j-1)*4+i, 14:23)))/max_norm,'b','LineWidth',2);    
       % plot(x_time, (aa((j-1)*5+i, 101:end)),'b','LineWidth',2);
        plot([aa((j-1)*4+i, 3)+3, aa((j-1)*4+i, 3)+3], [0,1], 'b', 'LineWidth',2);
        % correlation
        x_time_corr = 4:43;
        x_time_corr = x_time_corr + aa((j-1)*4+i,2) -1;
        plot(x_time_corr, veloc/max(veloc), 'r','LineWidth',2); 
        % start and end mark
        plot( x_start, y_marker, 'g-','LineWidth',1.5);
        plot( x_stop,  y_marker, 'g-','LineWidth',1.5);
        plot( [23.5,23.5],  y_marker, 'g-','LineWidth',1.5); % the middle
        plot( [29.5,29.5],  y_marker, 'g-','LineWidth',1.5); % the 300 ms delay
%         % set the same scale for all plot
         xlim([0,48]);
         ylim([0,1.05]);
    end 
%     axes('position',[0 0 1 1]); 
%     xlim([-50,50]);
%     ylim([-50,50]);
%     text(-25,45, 'vestibular'); 
%     text(0,45, 'visual'); 
%     text(25,45, 'combined'); 
%     text(-45, 45, FILE);
end
axis off;

%dlmwrite('alloutput.txt', peak_t');