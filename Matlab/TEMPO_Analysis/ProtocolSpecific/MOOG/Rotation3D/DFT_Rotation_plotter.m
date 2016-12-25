% This m file plots various quantities obtained by running
% Frequency_analysis, Rotation_frequency_analysis or frequency_analysis_2d.
% written by Narayan Ganesan (2007)
% REMEMBER to change the textread filenames to your specific files.(some notes - ab)

clear all;
str1 =  ' %s%f%f%f%f%f%f'; %6 values for direction specific (26/cell for the rotation and direction tuning protocols)
[name  dft_ratio, dft_p, dft_p_vel, dft_p_acc, phase, delay] = textread('Rotation_Frequency_26trajectories.dat', str1);
str1 =  ' %s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'; %22 values for cell specific
[name2 P_anova, P_anova_vel, P_anova_acc, DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, corr_vel_mfr, p_vel_mfr, corr_acc_mfr, p_acc_mfr, corr_acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc] = textread('Rotation_Frequency_cell.dat', str1);

vel_comp = -dft_ratio.*cos(phase);
acc_comp = dft_ratio.*sin(phase);

k=1; %index of cell_index

cell_index{k}(1) = 1;

for i=1:length(name)
    found=0;
    for i1 = 1:k
        if strcmp(name(cell_index{i1}(1)), name(i))
            if (i ~= 1)
                cell_index{i1}(length(cell_index{i1})+1) = i;
            end
            found=1;
        end
    end
    if (found==0)
        k=k+1;
        cell_index{k}(1) = i;
    end
end

% for i=1:length(cell_index)
%     dft_max_index(i) = find(abs(dft_ratio(cell_index{i}))==max(dft_ratio(cell_index{i})));
%     dft_max_names(i) = cell_index{i}(dft_max_index(i));
%     
%     vel_max_index(i) = find(abs(vel_comp(cell_index{i}))==max(vel_comp(cell_index{i})));
%     vel_max_names(i) = cell_index{i}(vel_max_index(i));
%     
%     acc_max_index(i) = find(abs(acc_comp(cell_index{i}))==max(acc_comp(cell_index{i})));
%     acc_max_names(i) = cell_index{i}(acc_max_index(i));
% end

%classify into significant velocity/accleration/both/none groups.
% and plot the DFTR vs. the no. of directions

total_dirs = length(dft_p);

for i=1:length(cell_index)
    sig_dirs(i) = length(find(strcmp(name,name(cell_index{i}(1))) & dft_p<0.05)); %number of significant directions for each cell
end

total_sig_dirs = sum(sig_dirs);

select_vel = find(P_anova_vel<0.05 & P_anova_acc>=0.05);
sig_vel = sum(sig_dirs(select_vel));
percent_vel = sig_vel/total_sig_dirs*100

select_acc = find(P_anova_vel>=0.05 & P_anova_acc<0.05);
sig_acc = sum(sig_dirs(select_acc));
percent_acc = sig_acc/total_sig_dirs*100

select_both = find(P_anova_vel<0.05 & P_anova_acc<0.05);
sig_both = sum(sig_dirs(select_both));
percent_both = sig_both/total_sig_dirs*100

select_none = find(P_anova_vel>=0.05 & P_anova_acc>=0.05);
sig_none = sum(sig_dirs(select_none));
percent_none = sig_none/total_sig_dirs*100

% figure
% hold on;
% dftr_binsize = 1%dftr binsize...
% edges = 1:dftr_binsize:10;%edges for histc computation
% 
% %plot significant acceleration or velocity, both or none
% sig_none_hist = histc(sig_none,edges); %compute histogram data
% bar(edges,sig_none_hist,0.2,'w');
% %plot significant velocity only
% sig_vel_hist = histc(sig_vel,edges);
% bar(edges+dftr_binsize/4,sig_vel_hist,0.2,'r'); %introduce offset to seperate the bars. offset = dftr_binsize/4
% %plot significant acceleration only
% sig_acc_hist = histc(sig_acc,edges);
% bar(edges+2*dftr_binsize/4,sig_acc_hist,0.2,'g');%offset = 2*dir_binsize/4
% %plot significant both
% sig_both_hist = histc(sig_both,edges);
% bar(edges+3*dftr_binsize/4,sig_both_hist,0.2,'k');%offset = 3*dir_binsize/4
% maxy = max([sig_none_hist sig_vel_hist sig_acc_hist sig_both_hist]);
% axis([0 10 0 maxy(1)+10]);
% title('DFTR vs no. of directions');
% %legend(num2str(percent_acc,3)+'%');
% legend(strcat('Neither(',num2str(percent_none,3),'%)'),strcat('Velocity(',num2str(percent_vel,3),'%)'),...
% strcat('Acceleration Only(',num2str(percent_acc,3),'%)'),strcat('Vel+Acc(',num2str(percent_both,3),'%)'));

%plot phase vs. the percentage of directions
figure;
edges=0:30:360;
temp_phase = phase*180/pi;
select = find(temp_phase<0);
temp_phase(select) = temp_phase(select) + 360; %wrap phase around from 0 to 360.
totalnoofdir = sum(histc(temp_phase, edges));
perctgdir = histc(temp_phase, edges)/totalnoofdir*100;
bar(edges, perctgdir,'w');
title('Phase vs. percentage of directions');
%set(H,'xtick',[0 60 120 180 240 300 360]);

%plot the scatter plot of VS Velocity, Acceleration, DDI velocty and
%Acceleration

%plot of no. of significant directions per cell..

unique_name = unique(name);
for i=1:length(unique_name)
    select_vel = find(strcmp(name, unique_name(i)) & dft_p_vel<0.05 & dft_p_acc>0.05);
    sig_dir_vel(i) = length(select_vel);%no. of significant directions for velocity
    select_both = find(strcmp(name, unique_name(i)) & dft_p_vel<0.05 & dft_p_acc<0.05);
    sig_dir_both(i) = length(select_both);
    select_acc = find(strcmp(name, unique_name(i)) & dft_p_vel>0.05 & dft_p_acc<0.05);
    sig_dir_acc(i) = length(select_acc);
    select_none = find(strcmp(name, unique_name(i)) & dft_p_vel>0.05 & dft_p_acc>0.05);
    sig_dir_none(i) = length(select_none);
end
figure
hold on;
dir_binsize = 2.5%direction binsize.. i.e no. of directions per bin.
edges = 1:dir_binsize:28;%bins for no. of directions

%plot no significant acceleration or velocity
bar(edges,histc(sig_dir_none,edges),0.2,'w');
%plot significant velocity only
bar(edges+dir_binsize/4,histc(sig_dir_vel,edges),0.2,'r'); %introduce offset to seperate the different color coded bars. offset = dir_binsize/4
%plot significant acceleration only
bar(edges+2*dir_binsize/4,histc(sig_dir_acc,edges),0.2,'g');%offset = 2*dir_binsize/4
%plot significant both
bar(edges+3*dir_binsize/4,histc(sig_dir_both,edges),0.2,'k');%offset = 3*dir_binsize/4
title('DFTR vs no. of directions');
%legend(num2str(percent_acc,3)+'%');
legend(strcat('Neither(',num2str(percent_none,3),'%)'),strcat('Velocity(',num2str(percent_vel,3),'%)'),...
strcat('Acceleration Only(',num2str(percent_acc,3),'%)'),strcat('Vel+Acc(',num2str(percent_both,3),'%)'));


figure
%Plots of Correlation coefficients between velocity and mfr
width=0.5;%width of bars
select = find(p_vel_mfr<0.05);
sig_vel_mfr = corr_vel_mfr(select);
select = find(p_vel_mfr>0.05);
nonsig_vel_mfr = corr_vel_mfr(select);
data1 = histc(sig_vel_mfr,-1.0:0.1:1.0);
data2 = histc(nonsig_vel_mfr,-1.0:0.1:1.0);
maxy=max([max(data1) max(data2)]);%find maximum y for plots
H = subplot(3,1,1);
bar(-1.0:0.1:1.0,data1,width,'k');
hold on;
bar(-1.0:0.1:1.0,data2,width,'w');
axis([-1 1 0 maxy(1)+1]);%+1 for buffer and need maxy(1) as dimension of max array might be more than 1
yticksteps = max([1,round(maxy(1)/4)]);%set y tick intervals
set(H,'ytick',0:yticksteps:maxy(1)+1);
set(H,'xtick',[-1 -0.5 0 0.5 1]);
set(H,'Box','Off');
set(H,'PlotBoxAspectRatio',[2,1,1]);
title('Velocity vs. MFR');
legend('p<0.05','p>0.05','location','best');
%set(H,'FontSize',30);

%Plots of Correlation coefficients between acceleration and mfr

select = find(p_acc_mfr<0.05);
sig_acc_mfr = corr_acc_mfr(select);
select = find(p_acc_mfr>0.05);
nonsig_acc_mfr = corr_acc_mfr(select);
data1 = histc(sig_acc_mfr,-1.0:0.1:1.0);
data2 = histc(nonsig_acc_mfr,-1.0:0.1:1.0);
maxy=max([max(data1) max(data2)]);
H = subplot(3,1,2);
bar(-1.0:0.1:1.0,data1,width,'k');
hold on;
bar(-1.0:0.1:1.0,data2,width,'w');
axis([-1 1 0 maxy(1)+1]);%+1 for buffer and need maxy(1) as dimension of max array might be more than 1
yticksteps = max([1,round(maxy(1)/4)]);%set y tick intervals
set(H,'ytick',0:yticksteps:maxy(1)+1);
set(H,'xtick',[-1 -0.5 0 0.5 1]);
set(H,'Box','Off');
set(H,'PlotBoxAspectRatio',[2,1,1]);
title('Acceleration vs. MFR');
legend('p<0.05','p>0.05','location','best');
%set(H,'FontSize',30);

%Plots of Correlation coefficients between acceleration and velocity

select = find(p_acc_vel<0.05);
sig_acc_vel = corr_acc_vel(select);
select = find(p_acc_vel>0.05);
nonsig_acc_vel = corr_acc_vel(select);
data1 = histc(sig_acc_vel,-1.0:0.1:1.0);
data2 = histc(nonsig_acc_vel,-1.0:0.1:1.0);
maxy=max([max(data1) max(data2)]);%find maximum y for plots
H = subplot(3,1,3);
bar(-1.0:0.1:1.0,data1,width,'k');
hold on;
bar(-1.0:0.1:1.0,data2,width,'w');
axis([-1 1 0 maxy(1)+1]);%+1 for buffer and need maxy(1) as dimension of max array might be more than 1
yticksteps = max([1,round(maxy(1)/4)]);%set y tick intervals
set(H,'ytick',0:yticksteps:maxy(1)+1);
set(H,'xtick',[-1 -0.5 0 0.5 1]);
set(H,'Box','Off');
set(H,'PlotBoxAspectRatio',[2,1,1]);
title('Acceleration vs. Velocity');
legend('p<0.05','p>0.05','location','best');
%set(H,'FontSize',30);
