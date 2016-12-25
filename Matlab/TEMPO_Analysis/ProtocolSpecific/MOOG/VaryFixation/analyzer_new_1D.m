%analyze data stored in fit_data which contains 16 variables.
%When loading fit_data, or the data for the visual condition, change the
%criterion so that each cell has 5 directions at least (except for some
%figures)
%separate the figures into the 5 groups A1, A2L, A2R, A3 and V as opposed
%to by significance.

clear all;
% filepath = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\VaryFixation\';

% for n = 1:3
n = 1;

clear FILE azimuth gaze p_val VAF_com VAF_vel DC_com b_com tau_com sigma_com a_com DC_vel K_vel tau_vel sigma_vel DC_unc K_unc tau_unc sigma_unc DFT_ratio;

if n == 1
    filename = 'PSTH_fitdata_ves.mat';
elseif n == 2
    filename = 'PSTH_fitdata_vis.mat';
else
    filename = 'PSTH_fitdata_comb.mat';
end

[FILE azimuth gaze p_val VAF_com VAF_vel DC_com b_com tau_com sigma_com a_com DC_vel K_vel tau_vel sigma_vel DC_unc K_unc tau_unc sigma_unc DFT_ratio] = textread([filepath filename], ... 
'%s %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', '\t', 'headerlines', 1);

% OR

load(filename);

%--------------------------------------------------------------------------
% Find the distribution of the VAFs

k = 1;
for i = .05:.05:1
	VAF_com_distribution{n}(k) = length(find( (VAF_com < i) & (VAF_com >= i - 0.05) ));
	VAF_vel_distribution{n}(k) = length(find( (VAF_vel < i) & (VAF_vel >= i - 0.05) ));
	k = k+1;
end

figure(10*n)
bar(0:.05:.95, VAF_com_distribution{n})
title('Combined Model')
xlabel('VAF')
ylabel('Frequency')

figure(10*n+1)
bar(0:.05:.95, VAF_vel_distribution{n})
title('Velocity Only Model')
xlabel('VAF')
ylabel('Frequency')


%--------------------------------------------------------------------------
%Find the distribution of the P Values

%Bin p values into bins of 0.01. Then the first bin will show the number of
%directions that will be significant.(total no. of directions is 2088)
k=1;
for i = 0:.01:1
    p_distribution{n}(k) = length(find( (p_val < i) & (p_val >= i - .01) ));
    k = k+1;
end

figure(10*n+3)
bar(0:.01:1, p_distribution{n})
title('Distribution of P Values')
xlabel('P value')
ylabel('Frequency')


%--------------------------------------------------------------------------
%Find the ratio of the velocity component to the absolute value of the acceleration component. 
%Also find the ratio of the DFT as a function of direction.

unique_azimuth = unique(azimuth);
unique_condition_num = unique(gaze);
b_to_a{n} = b_com./a_com;

k=1;
for i = -20:.1:20
    bbb(k) = length(find(b_to_a{n} >= i & b_to_a{n} < i+.1));
    k = k+1;
end


for j=1:length(unique_condition_num)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth==unique_azimuth(i)) & (gaze==unique_condition_num(j)) );
        b_to_a_ratio{n,i,j} = b_to_a{n}(select);
        tuned_dir{n,i,j} = DFT_ratio(select);
    end
end
%--------------------------------------------------------------------------

%Find the mean for each direction

for j=1:length(unique_condition_num)
    for i=1: length(unique_azimuth)
        %take the median to ignore outliers. outliers occur when a is
        %almost equal to zero.
        median_ba_ratio{n,i,j} = median(b_to_a_ratio{n,i,j});
        mean_tuned_dir{n,i,j} = mean(tuned_dir{n,i,j});
    end
end

%first we throw out directions which have VAF's of less than 60% in the
%combined model.

directions = find(VAF_com >= 0.6);
good_dir = FILE(directions);

%Define the max direction of each cell as being the direction with the
%highest DFT ratio (cell shows the greatest tuning). We only will select
%cells which had at least 5 directions which were selected to be fitted and
%which passed the VAF criterion.


unique_name = unique(good_dir);
match = zeros(1,length(unique_name));

%when you use unique, it sorts the names so they are no longer in order.
%This could cause problems when indexing later. To keep the same order, we
%must reorder the unique_name vector

unique_index = zeros(1,length(pref_cell_name));
for i = 1:length(pref_cell_name)
    for j = 1:length(unique_name)
        if (length(unique_name{j}) == length(pref_cell_name{i})) & (sum(unique_name{j} == pref_cell_name{i}) == length(unique_name{j}))
            unique_index(i) = 1;
        end
    end
end
unique_index = find(unique_index == 1);
unique_name = pref_cell_name(unique_index);

for i = 1:length(unique_name)
    for j = 1:length(good_dir)
        i
        if (length(unique_name{i}) == length(good_dir{j})) & (sum(unique_name{i} == good_dir{j}) == length(unique_name{i}))
            match(i) = match(i) + 1;
        end
        if match(i) >=1
            matches(i) = 1;
        else
            matches(i) = 0;
        end
    end
end
matches = find(matches == 1);


%these cells have at least 5 directions which were selected by the DFT
%method (ratio was greater than 3).

good_cells = unique_name(matches);

new_DFT_ratio = DFT_ratio(directions);
good_azimuth = azimuth(directions);
good_gaze = gaze(directions);
good_tau_com = tau_com(directions);
good_sig_com = sigma_com(directions);


%The max direction will be the direction which has the highest DFT ratio.
for i = 1:length(good_cells)
    max_ratio = 0;
    for j = 1:length(good_dir)
        if (length(good_cells{i}) == length(good_dir{j})) & (sum(good_cells{i} == good_dir{j}) == length(good_cells{i})) & (max_ratio < new_DFT_ratio(j))
            max_ratio = new_DFT_ratio(j);
            max_az = good_azimuth(j);
            max_gaze = good_gaze(j);
            max_index(i) = j;
        end
    end
    %max_dir contains the max [azimuth, Elevation] for each of the good
    %cells. Its index corresponds to that of good_cells such that
    %the max direction of good_cells{1} is max_dir{1}.
    
    max_dir{i} = [max_az, max_gaze];
    max_azimuth(i) = max_az;
    max_gaze(i) = max_gaze;
end
tau_max_com = good_tau_com(max_index);


% %find angle between preferred direction and other directions
% Azims = [0 45 90 135 180 225 270 315]*pi/180;
% Elevs = [-90 -45 0 45 90]*pi/180;

%b/a ratio of each of the good directions
ba_vals = b_to_a{n}(directions);

angle =[];
diff_angle{length(good_cells)} = [];
ba_rat{length(good_cells)} = [];

direc = zeros(0, length(good_dir));
count =1;
 
% for i = 1:length(good_cells)
%     for j = 1:length(good_dir)
%         if (length(good_cells{i}) == length(good_dir{j})) & (sum(good_cells{i} == good_dir{j}) == length(good_cells{i}))
%             temp_max_az = Azims(max_azimuth(i));
%             temp_max_el = Elevs(max_elevation(i));
%             temp_azimuth = Azims(good_azimuth(j));
%             temp_elevation = Elevs(good_elevation(j));
%             direc(count) = j;
%             count = count+1;
%             % calculate difference angles using the formula:
%             % cos(x) = sin(Ea)sin(Eb) + cos(Eb)sin(Ab)cos(Ea)sin(Aa) + cos(Eb)cos(Ab)cos(Ea)cos(Aa)
%             % where Ea and Eb is the elevation of vectors a and b, and Aa and Ab =
% %             azimuth of a and b and x is the difference angle
%             angle = round(acos(sin(temp_max_el)*sin(temp_elevation)+cos(temp_elevation)*sin(temp_azimuth)*cos(temp_max_el)*sin(temp_max_az)+cos(temp_elevation)*cos(temp_azimuth)*cos(temp_max_el)*cos(temp_max_az))*180/pi);
%             
%             diff_angle{i} = [diff_angle{i} angle];
%             ba_rat{i} = [ba_rat{i} ba_vals(j)];            
%         end
%     end
% end
% 
% 
% 
% for i = 1:length(diff_angle)
%     unique_angle{i} = unique(diff_angle{i});
%     
%     for j = 1:length(unique_angle{i})
%         find_index = find(diff_angle{i} == unique_angle{i}(j));
%         good_ba_ratio{i}(j) = median(ba_rat{i}(find_index));
%     end
% end


%scatter plot that Dora wanted with just the max directions (for this to
%work, remember to select all cells that have AT LEAST 1 direction fitted.
%Not AT LEAST 5 like we have been doing for everything else). The figure
%from this is saved as scatterplot(max directions).fig so you dont have to
%run this part at all.

%Scatter plot

good_VAF_com = VAF_com(directions);
good_VAF_vel = VAF_vel(directions);

%good_com_dir are the VAFs for the combined model directions that
%passed both the VAF and DFT criteria

good_veloc_dir = good_VAF_vel(direc);
good_com_dir = good_VAF_com(direc);
good_p1 = p_val(directions);
good_p = good_p1(direc);
significant = find(good_p <= 0.05);
% signif2 = find(p <= 0.05);

max_veloc_dir = good_VAF_vel(max_index);
max_com_dir = good_VAF_com(max_index);
max_p = good_p1(max_index);
signifi = find(max_p <= 0.05);
not_signifi = find(max_p > 0.05);
%find good cells which are significant with p<.05
significant_cells = find(max_p<=0.05);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Define different populations as being:
%A1 = significant cells with 1< free tau <2
%A2L = significant cells with 2<= free tau < 2.25
%A2R = significant cells with 2.25<= free tau <2.5
%A3 = significant cells with 2.5<= free tau <=3
%V = non-signifcicant cells

%find the tau's for the max_directions
free_tau = tau_unc(directions);
free_tau = free_tau(max_index);

k1=1;k2=1;k3=1;k4=1;k6=1;
for i = 1:length(max_index)
    if (free_tau(i) >= 0 & free_tau(i) < 2) & max_p(i) <= 0.05
        A1(k1) = i;
        k1=k1+1;
    elseif (free_tau(i) >= 2 & free_tau(i) < 2.25) & max_p(i) <= 0.05
        A2L(k2) = i;
        k2=k2+1;
    elseif (free_tau(i) >= 2.25 & free_tau(i) < 2.5) & max_p(i) <= 0.05
        A2R(k6) = i;
        k6=k6+1;
    elseif (free_tau(i) >= 2.5 & free_tau(i) <= 3) & max_p(i) <= 0.05
        A3(k3) = i;
        k3=k3+1;
    elseif free_tau(i) > 1 & max_p(i) > 0.05
        V(k4) = i;
        k4=k4+1;
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% figure
% dum_vel2 = max_com_dir;
% plot(dum_vel2, max_com_dir, 'y', 'linewidth', 2)
% hold on
% plot(max_veloc_dir(V), max_com_dir(V), 'g^', 'markerfacecolor', 'g')
% plot(max_veloc_dir(A1), max_com_dir(A1), 'bs', 'markerfacecolor', 'b') 
% xlabel('VAF velocity')
% ylabel('VAF combined')
% plot(max_veloc_dir(A2L), max_com_dir(A2L), 'ro', 'markerfacecolor', 'r')
% plot(max_veloc_dir(A2R), max_com_dir(A2R), 'ko', 'markerfacecolor', 'k')
% plot(max_veloc_dir(A3), max_com_dir(A3), 'ms', 'markerfacecolor', 'm')
% 
% title('Scatter plot of VAF from combined and velocity only models for max directions')
% legend('Unity', 'V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% 
% %--------------------------------------------------------------------------
% 
% %Find the matches of the good cells with Yong's cells; the cells which
% %passed the VAF, 5 directions and DFT criteria
% 
% max_match = zeros(1,length(yong_cell_name));
% for i = 1:length(yong_cell_name)
%     for j = 1:length(good_cells)
%         if (length(yong_cell_name{i}) == length(good_cells{j})) & (sum(yong_cell_name{i} == good_cells{j}) == length(yong_cell_name{i}))
%             max_match(i) = 1;
%         end
%     end
% end
% max_match = find(max_match == 1);
% 
% yong_max_matched_cells = yong_cell_name(max_match);
% 
% %good cells matched with the pref dirs from the DFT method and yongs
% %method. The pref direction is different from the max direction. The max
% %direction is the direction that was recorded from which had a maximal
% %response. The preferred direction is the direction calculated from the
% %vector sum method. 
% DFT_pref_matched = pref_cell_name(max_match);
% DFT_pref_az = pref_az(max_match);
% DFT_pref_el = pref_el(max_match);
% 
% yong_pref_matched = yong_pref_cell_name(max_match);
% yong_pref_az = yong_pref_az(max_match);
% yong_pref_el = yong_pref_el(max_match);
% 
% %max directions from yongs method
% yong_max_el = yong_max_el(max_match);
% yong_max_az = yong_max_az(max_match);
% %--------------------------------------------------------------------------
% 
% 
% %plot b/a ratio vs tau for max directions, separate by population
% badum = b_to_a(directions);
% b_a_max = badum(max_index);
% 
% figure
% plot(b_a_max, b_a_max, 'y', 'linewidth', 2)
% hold on
% plot(tau_max_com(V), b_a_max(V), 'g^', 'markerfacecolor', 'g')
% plot(tau_max_com(A1), b_a_max(A1),'bs', 'markerfacecolor', 'b')
% plot(tau_max_com(A2L), b_a_max(A2L), 'ro', 'markerfacecolor', 'r')
% plot(tau_max_com(A2R), b_a_max(A2R), 'ko', 'markerfacecolor', 'k')
% plot(tau_max_com(A3), b_a_max(A3), 'ms', 'markerfacecolor', 'm')
% 
% axis([1.5 3 -20 20])
% xlabel('Tau')
% ylabel('b/a ratio')
% title('b/a ratio vs. tau(combined) for max directions')
% legend('Unity' , 'V', 'A1', 'A2L', 'A2R', 'A3')
% 
% 
% 
% %plot of figure from yong's paper of pref azimuth vs pref elevation
% 
% figure
% 
% plot(DFT_pref_az(V), DFT_pref_el(V), 'g^', 'markerfacecolor', 'g')
% hold on
% plot(DFT_pref_az(A1), DFT_pref_el(A1),'bs', 'markerfacecolor', 'b')
% plot(DFT_pref_az(A2L), DFT_pref_el(A2L), 'ro', 'markerfacecolor', 'r')
% plot(DFT_pref_az(A2R), DFT_pref_el(A2R), 'ko', 'markerfacecolor', 'k')
% plot(DFT_pref_az(A3), DFT_pref_el(A3), 'ms', 'markerfacecolor', 'm')
% 
% title('Plot of Preferred Azimuth vs Preferred Elevation')
% xlabel('Preferred Azimuth')
% ylabel('Preferred Elevation')
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% %--------------------------------------------------------------------------
% %plot the mean b/a ratio against preferred azimuth and elevation from the
% %combined method
% for i = 1:length(ba_rat)
%     temp_ba_rat = ba_rat{i}(ba_rat{i} < 25);
% mean_ba(i) = mean(temp_ba_rat);
% end
% 
% %wrap
% 
% DFT_pref_az(DFT_pref_az >= 270) = DFT_pref_az(DFT_pref_az >= 270) - 360;        
%  
% figure
% plot(DFT_pref_az(V), mean_ba(V), 'g^', 'markerfacecolor', 'g')
% hold on
% plot(DFT_pref_az(A1), mean_ba(A1), 'bs', 'markerfacecolor', 'b')
% plot(DFT_pref_az(A2L), mean_ba(A2L), 'ro', 'markerfacecolor', 'r')
% plot(DFT_pref_az(A2R), mean_ba(A2R), 'ko', 'markerfacecolor', 'k')
% plot(DFT_pref_az(A3), mean_ba(A3), 'ms', 'markerfacecolor', 'm')
% 
% xlabel('Preferred Azimuth')
% ylabel('Mean b/a ratio')
% title('Mean b/a ratio vs. preferred azimuth')
% axis([-90 270 -10 10])
% set(gca,'xtick', [-90 -45 0 45 90 135 180 225 270])
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% 
% figure
% plot(DFT_pref_el(V), mean_ba(V), 'g^', 'markerfacecolor', 'g')
% hold on
% plot(DFT_pref_el(A1), mean_ba(A1), 'bs', 'markerfacecolor', 'b')
% plot(DFT_pref_el(A2L), mean_ba(A2L), 'ro', 'markerfacecolor', 'r')
% plot(DFT_pref_el(A2R), mean_ba(A2R), 'ko', 'markerfacecolor', 'k')
% plot(DFT_pref_el(A3), mean_ba(A3), 'ms', 'markerfacecolor', 'm')
% 
% 
% xlabel('Preferred Elevation')
% ylabel('Mean b/a ratio')
% title('Mean b/a ratio vs. preferred elevation')
% axis([-90 90 -10 10])
% set(gca, 'xtick', [-90 -45 0 45 90])
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% DFT_pref_az = pref_az(max_match);
% %--------------------------------------------------------------------------       
% %make a scatter plot of sigma and tau for the max directions only
% tau_max_dir = good_tau_com(max_index);
% sig_max_dir = good_sig_com(max_index);
% dum_s = tau_max_dir;
% 
% figure
% plot(tau_max_dir(V), sig_max_dir(V), 'g^', 'markerfacecolor', 'g')
% hold on
% plot(tau_max_dir(A1), sig_max_dir(A1), 'bs', 'markerfacecolor', 'b')
% plot(tau_max_dir(A2L), sig_max_dir(A2L), 'ro', 'markerfacecolor', 'r')
% plot(tau_max_dir(A2R), sig_max_dir(A2R), 'ko', 'markerfacecolor', 'k')
% plot(tau_max_dir(A3), sig_max_dir(A3), 'ms', 'markerfacecolor', 'm')
% 
% title('Sigma vs Tau for the max directions')
% xlabel('tau')
% ylabel('sigma')
% axis([0 3 0 1])
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% %--------------------------------------------------------------------------
% 
% %distribution of free tau
% good_tau_vel = tau_velocity(directions);
% tau_max_vel = good_tau_vel(max_index);
% ss = find(p<= 0.05);
% nn = find(p> 0.05);
% max_sig = find(max_p <= 0.05);
% max_not = find(max_p > 0.05);
% %for cells fitted with only velocity model and tau unrestrained
% k=1;
% for i = 0:.05:3
%     V_tau_dist(k) = length(find((tau_max_vel(V) >= i) & (tau_max_vel(V) < i + 0.05)));
%     A1_tau_dist(k) = length(find((tau_max_vel(A1) >= i) & (tau_max_vel(A1) < i + 0.05)));
%     A2L_tau_dist(k) = length(find((tau_max_vel(A2L) >= i) & (tau_max_vel(A2L) < i + 0.05)));
%     A2R_tau_dist(k) = length(find((tau_max_vel(A2R) >= i) & (tau_max_vel(A2R) < i + 0.05)));
%     A3_tau_dist(k) = length(find((tau_max_vel(A3) >= i) & (tau_max_vel(A3) < i + 0.05)));   
%     k = k+1;
% end
% 
% figure
% bar(0:.05:3, V_tau_dist,.5,'g')
% hold on
% bar(0.025:.05:3.025, A1_tau_dist,.5,'b')
% bar(0.025:.05:3.025, A2L_tau_dist,.5,'r')
% bar(0.025:.05:3.025, A2R_tau_dist,.5,'k')
% bar(0:.05:3, A3_tau_dist,.5,'m')
% xlabel('Tau')
% ylabel('Frequency')
% title('Distribution of Tau for max directions')
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% 
% % % orient landscape
% % % print -dwinc
% % % close
% % 
% 
% %--------------------------------------------------------------------------
% % plot scatter plot of DFT preferred direction vs Yongs preferred direction
% 
% % wrap data
% for i = 1:length(DFT_pref_az)
%     if DFT_pref_az(i) - yong_pref_az(i) > 180 & DFT_pref_az(i) > yong_pref_az(i)
%         DFT_pref_az(i) = DFT_pref_az(i) - 360;
%     end
%     if yong_pref_az(i) - DFT_pref_az(i) > 180 & yong_pref_az(i) > DFT_pref_az(i)
%         yong_pref_az(i) = yong_pref_az(i) - 360;
%     end
% end
%     
% figure
% plot(yong_pref_az, yong_pref_az, 'y', 'linewidth', 2)
% hold on
% plot(yong_pref_az(V), DFT_pref_az(V), 'g^', 'markerfacecolor', 'g')
% plot(yong_pref_az(A1), DFT_pref_az(A1), 'bs', 'markerfacecolor', 'b')
% plot(yong_pref_az(A2L), DFT_pref_az(A2L), 'ro', 'markerfacecolor', 'r')
% plot(yong_pref_az(A2R), DFT_pref_az(A2R), 'ko', 'markerfacecolor', 'k')
% plot(yong_pref_az(A3), DFT_pref_az(A3), 'ms', 'markerfacecolor', 'm')
% title('Scatter plot of Preferred Azimuths from DFT and Mean Firing rate methods')
% xlabel('Pref Azimuth (Mean Firing rate method)')
% ylabel('Pref Azimuth (DFT method)')
% axis([-135 360 -135 360])
% set(gca, 'xtick', -180:45:360)
% set(gca, 'ytick', -180:45:360)
% legend('Unity' , 'V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% 
% figure
% plot(yong_pref_el, yong_pref_el, 'y', 'linewidth', 2)
% hold on
% plot(yong_pref_el(V), DFT_pref_el(V), 'g^', 'markerfacecolor', 'g')
% plot(yong_pref_el(A1), DFT_pref_el(A1), 'bs', 'markerfacecolor', 'b')
% plot(yong_pref_el(A2L), DFT_pref_el(A2L), 'ro', 'markerfacecolor', 'r')
% plot(yong_pref_el(A2R), DFT_pref_el(A2R), 'ko', 'markerfacecolor', 'k')
% plot(yong_pref_el(A3), DFT_pref_el(A3), 'ms', 'markerfacecolor', 'm')
% title('Scatter plot of Preferred Elevations from DFT and Mean Firing rate methods')
% xlabel('Pref Elevation (Mean Firing rate method)')
% ylabel('Pref Elevation (DFT method)')
% axis([-90 90 -90 90])
% set(gca, 'xtick', -90:45:90)
% set(gca, 'ytick', -90:45:90)
% legend('Unity' , 'V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% %--------------------------------------------------------------------------
% %find the angle between mean firing rate pref dir and dft pref dir
% 
% for i = 1:length(max_match)
%     % calculate difference angles using the formula:
%     % cos(x) = sin(Ea)sin(Eb) + cos(Eb)sin(Ab)cos(Ea)sin(Aa) + cos(Eb)cos(Ab)cos(Ea)cos(Aa)
%     % where Ea and Eb is the elevation of vectors a and b, and Aa and Ab =
%     % azimuth of a and b and x is the difference angle
%     Aa = yong_pref_az(i)*(pi/180);
%     Ab = DFT_pref_az(i)*(pi/180);
%     Ea = yong_pref_el(i)*(pi/180);
%     Eb = DFT_pref_el(i)*(pi/180);
%     direction_angle(i) = (acos(sin(Ea)*sin(Eb) + cos(Eb)*sin(Ab)*cos(Ea)*sin(Aa) + cos(Eb)*cos(Ab)*cos(Ea)*cos(Aa))*180/pi);
% end
% sig_direction_angle = direction_angle(signifi);
% not_sig_direction_angle = direction_angle(not_signifi);
% %--------------------------------------------------------------------------
% %plot direction_angle as a function of b/a ratio
% 
% max_direction_ba_ratio = ba_vals(max_index);
% sig_max_dir_ba = max_direction_ba_ratio(signifi);
% not_sig_max_dir_ba = max_direction_ba_ratio(not_signifi);
% 
% figure
% plot(direction_angle(V), max_direction_ba_ratio(V), 'g^', 'markerfacecolor', 'g')
% hold on
% plot(direction_angle(A1), max_direction_ba_ratio(A1), 'bs', 'markerfacecolor', 'b')
% plot(direction_angle(A2L), max_direction_ba_ratio(A2L), 'ro', 'markerfacecolor', 'r')
% plot(direction_angle(A2R), max_direction_ba_ratio(A2R), 'ko', 'markerfacecolor', 'k')
% plot(direction_angle(A3), max_direction_ba_ratio(A3), 'ms', 'markerfacecolor', 'm')
% title('Plot of difference angle vs. b/a ratio')
% xlabel('Difference Angle')
% ylabel('b/a ratio')
% axis([0 180 -20 20])
% legend('V', 'A1', 'A2L', 'A2R', 'A3')
% % orient landscape
% % print -dwinc
% % close
% 
% %_-------------------------------------------------------------------------
% 
% %get HTI figures
% load HTI_data_new
% k=1;
% %serparate by significance of HTI's
% for i = 0:.05:.6
%     V_vis(k) = length(find(HTI_vis(V)>= i & HTI_vis(V) < i+.05));
%     A1_vis(k) = length(find(HTI_vis(A1)>= i & HTI_vis(A1) < i+.05));
%     A2L_vis(k) = length(find(HTI_vis(A2L)>= i & HTI_vis(A2L) < i+.05));
%     A2R_vis(k) = length(find(HTI_vis(A2R)>= i & HTI_vis(A2R) < i+.05));
%     A3_vis(k) = length(find(HTI_vis(A3)>= i & HTI_vis(A3) < i+.05));
%     k=k+1;
% end
% 
% HTI_vis =HTI_vis(max_match);
% HTI_vest = HTI_vest(max_match);
% %seperate by significance of fit model.
% k =1;
% for i = 0:.05:.6
%     visual_sig_fit(k) = length(find(HTI_vis(signifi)>= i & HTI_vis(signifi) < i+.05));
%     visual_not_sig_fit(k) = length(find(HTI_vis(not_signifi) >=i & HTI_vis(not_signifi) <i+.05));
%     vest_sig_fit(k) = length(find(HTI_vest(signifi) >=i & HTI_vest(signifi) <i+.05));
%     vest_not_sig_fit(k) = length(find(HTI_vest(not_signifi) >=i & HTI_vest(not_signifi) <i+.05));
%     k=k+1;
% end
% 
% % figure
% % bar(0:.05:.6, binned_vis_sig,.5,'r')
% % hold on
% % bar(0.025:.05:.625,binned_vis_bad, .5) 
% % axis([0 .6 0 max([max(binned_vis_bad) max(binned_vis_sig)])])
% % xlabel('HTI')
% % ylabel('Number of Neurons')
% % title('Visual Condition')
% % legend('p<.05', 'p>.05')
% % orient landscape
% % print -dwinc
% % close
% 
% 
% % figure
% % bar(0:.05:.6, binned_vest_sig, .5,'r')
% % hold on
% % bar(0.025:.05:.625, binned_vest_bad, .5)
% % axis([0 .6 0 max([max(binned_vest_bad) max(binned_vest_sig)])])
% % xlabel('HTI')
% % ylabel('Number of Neurons')
% % title('Vestibular Condition')
% % legend('p<.05', 'p>.05')
% % orient landscape
% % print -dwinc
% % close
% 
% % figure
% % bar(0.0125:.05:.6125, visual_sig_fit,.5,'r')
% % hold on
% % bar(0.0375:.05:.6375,visual_not_sig_fit, .5) 
% % axis([0 .6 0 max([max(visual_not_sig_fit) max(visual_sig_fit)])])
% % xlabel('HTI distribution for max directions')
% % ylabel('Number of Neurons')
% % title('HTI distribution for max directions: Visual Condition')
% % legend('p<.05', 'p>.05')
% % orient landscape
% % print -dwinc
% % close
% 
% 
% % figure
% % bar(0.0125:.05:.6125, vest_sig_fit, .5,'r')
% % hold on
% % bar(0.0375:.05:.6375, vest_not_sig_fit, .5)
% % axis([0 .6 0 max([max(vest_not_sig_fit) max(vest_sig_fit)])])
% % xlabel('HTI')
% % ylabel('Number of Neurons')
% % title('HTI distribution for max directions: Vestibular Condition')
% % legend('p<.05', 'p>.05')
% % orient landscape
% % print -dwinc
% % close
% 
% 
% %--------------------------------------------------------------------------
% %to use this part of the code, first run analyzer.m with the first few
% %lines uncommented and the rest commented, to generate ba_vestibular and
% %good_vestibular cells. Remember to change the data that is being loaded
% %aswell
% 
% % -------------------------------------------------------------
% % ba_vestibular = ba_vals(max_index);
% % good_vestibular_cells = good_cells;
% % vestibular_azimuth = azimuth(directions);
% % vestibular_elevation = elevation(directions);
% % max_index_vestibular = max_index;
% % max_match_vestibular = max_match;
% % free_tau_vest = free_tau;
% % max_p_vest = max_p;
% % load HTI_data
% % index_vis = zeros(1,255);
% % index_vest = zeros(1,255);
% % index_vest(find(p_vest <= 0.05)) = 1;
% % index_vis(find(p_vis <=0.05)) = 1;
% % index = index_vest+index_vis;
% % index = find(index == 2);
% % vest_az = pref_az(index);
% % vest_el = pref_el(index);
% % direction_angle_vest = direction_angle;
% %------------------------------------------------------------- 
% 
% 
% 
% 
% % % %we need to find the cells that have been fitted in both the vestibular and
% % % %visual conditions in the max direction which is defined as the visual max
% % % %direction.
% % 
% % vis_vest_match = zeros(1, length(good_cells));
% % for i = 1:length(good_vestibular_cells)
% %     for j = 1:length(good_cells)
% %         if (length(good_vestibular_cells{i}) == length(good_cells{j})) & (sum(good_vestibular_cells{i} == good_cells{j}) == length(good_vestibular_cells{i}))
% %             vis_vest_match(j) = 1;
% %         end
% %     end
% % end
% % 
% % matched_cell_index = find(vis_vest_match == 1);
% % 
% % dummy_cells = good_cells(matched_cell_index);
% % for i = 1:length(dummy_cells)
% %     for j = 1:length(good_vestibular_cells)
% %         if (length(dummy_cells{i}) == length(good_vestibular_cells{j})) & (sum(dummy_cells{i} == good_vestibular_cells{j}) == length(dummy_cells{i}))
% %             vis_vest_match2(j) = 1;
% %         end
% %     end
% % end
% % matched_cell_index2 = find(vis_vest_match2 == 1);
% % 
% % free_tau = free_tau(matched_cell_index);
% % max_p = max_p(matched_cell_index);
% % k1=1;k2=1;k3=1;k4=1;k6=1;
% % for i = 1:length(matched_cell_index)
% %     if (free_tau(i) > 1 & free_tau(i) < 2) & max_p(i) <= 0.05
% %         A1_com(k1) = i;
% %         k1=k1+1;
% %     elseif (free_tau(i) >= 2 & free_tau(i) < 2.25) & max_p(i) <= 0.05
% %         A2L_com(k2) = i;
% %         k2=k2+1;
% %     elseif (free_tau(i) >= 2.25 & free_tau(i) < 2.5) & max_p(i) <= 0.05
% %         A2R_com(k6) = i;
% %         k6=k6+1;
% %     elseif (free_tau(i) >= 2.5 & free_tau(i) <= 3) & max_p(i) <= 0.05
% %         A3_com(k3) = i;
% %         k3=k3+1;
% %     elseif free_tau(i) > 1 & max_p(i) > 0.05
% %         V_com(k4) = i;
% %         k4=k4+1;
% %     end
% % end
% % 
% % ba_visual = ba_vals(max_index);
% % ba_visual = ba_visual(matched_cell_index);
% % 
% % ba_vestibular = ba_vestibular(matched_cell_index2);
% % % 
% % % 
% % % % figure
% % % % plot(ba_visual, ba_visual, 'y', 'linewidth', 2)
% % % % hold on
% % % % plot(ba_visual(A1_com), ba_vestibular(A1_com), 'ro', 'markerfacecolor', 'r')
% % % % plot(ba_visual(A2L_com), ba_vestibular(A2L_com), 'go', 'markerfacecolor', 'g')
% % % % plot(ba_visual(A2R_com), ba_vestibular(A2R_com), 'co', 'markerfacecolor', 'c')
% % % % plot(ba_visual(A3_com), ba_vestibular(A3_com), 'ko', 'markerfacecolor', 'k')
% % % % plot(ba_visual(V_com), ba_vestibular(V_com), 'bo', 'markerfacecolor', 'b')
% % % % plot(ba_visual(N_com), ba_vestibular(N_com), 'mo', 'markerfacecolor', 'm')
% % % % axis([-15 15 -15 15])
% % % % title('Scatter plot for b/a ratios for visual and vestibular conditions for all max directions')
% % % % xlabel('Visual b/a ratio')
% % % % ylabel('Vestibular b/a ratio')
% % % % legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V', 'N')
% % % %--------------------------------------------------------------------------
% % distribution of diff angle from vis/vest conditions
% % 
% % % visual_az = pref_az(index);
% % % visual_el = pref_el(index);
% % % for i = 1:length(visual_az)
% % %     Aa = visual_az(i)*(pi/180);
% % %     Ab = vest_az(i)*(pi/180);
% % %     Ea = visual_el(i)*(pi/180);
% % %     Eb = vest_el(i)*(pi/180);
% % %     dif_ang_vis_vest(i) = (acos(sin(Ea)*sin(Eb) + cos(Eb)*sin(Ab)*cos(Ea)*sin(Aa) + cos(Eb)*cos(Ab)*cos(Ea)*cos(Aa))*180/pi);
% % % end
% % % 
% % % % plot a distribution of the delta preferred direction
% % % k = 1;
% % % for i = 0:15:180
% % %     binned_angle(k) = length(find(dif_ang_vis_vest >= i & dif_ang_vis_vest < i+15));
% % %     k = k+1;
% % % end
% % % figure
% % % bar(0:15:180, binned_angle)
% % % title('Distribution of the difference between the Pref Directions from the visual and vestibular conditions')
% % % xlabel('Difference in angle between conditions')
% % % ylabel('Number')
% % % set(gca, 'xtick', 0:15:180)
% % %--------------------------------------------------------------------------
% % %scatter plot of diff angle from vis/vest conditions
% % 
% % direction_angle_vis = direction_angle(matched_cell_index);
% % direction_angle_vest = direction_angle_vest(matched_cell_index2);
% % figure
% % plot(direction_angle_vis, direction_angle_vis, 'y', 'linewidth', 2)
% % hold on
% % plot(direction_angle_vis(A1_com),direction_angle_vest(A1_com),'o', 'markerfacecolor', 'b')
% % plot(direction_angle_vis(A2L_com),direction_angle_vest(A2L_com),'ro', 'markerfacecolor', 'r')
% % plot(direction_angle_vis(A2R_com),direction_angle_vest(A2R_com),'go', 'markerfacecolor', 'g')
% % plot(direction_angle_vis(A3_com),direction_angle_vest(A3_com),'k^', 'markerfacecolor', 'k')
% % plot(direction_angle_vis(V_com),direction_angle_vest(V_com),'m^', 'markerfacecolor', 'm')
% % legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V')
% % title('Scatter plot of diff angle btwn mean firing rate and DFT method- visual vs. vest')
% % xlabel('Visual Difference Angle')
% % ylabel('Vestibular Difference Angle')
% 
% 
end