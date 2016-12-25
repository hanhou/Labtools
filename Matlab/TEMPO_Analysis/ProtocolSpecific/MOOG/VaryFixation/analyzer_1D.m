%analyze data stored in fit_data which contains 16 variables.
%When loading fit_data, or the data for the visual condition, change the
%criterion so that each cell has 5 directions at least (except for some
%figures)

load fit_data_chris

% load fit_data;

% load fit_data_vestibular

% Find the distribution of the VAFs

k = 1;
for i = .05:.05:1
VAF_com_distribution(k) = length(find( (VAF_com < i) & (VAF_com >= i - 0.05) ));
VAF_vel_distribution(k) = length(find( (VAF_vel < i) & (VAF_vel >= i - 0.05) ));
k = k+1;
end
% figure
% bar(0:.05:.95, VAF_com_distribution)
% title('Combined Model')
% xlabel('VAF')
% ylabel('Frequeny')
% 
% figure
% bar(0:.05:.95, VAF_vel_distribution)
% title('Velocity Only Model')
% xlabel('VAF')
% ylabel('Frequeny')
%--------------------------------------------------------------------------
%Find the distribution of the P Values

for i = 1:length(p_val)
    p(i) = str2num(p_val{i});
end

%Bin p values into bins of 0.01. Then the first bin will show the number of
%directions that will be significant.(total no. of directions is 2088)
k=1;
for i = 0:.01:1
    p_distribution(k) = length(find( (p < i) & (p >= i - .01) ));
    k = k+1;
end

% figure
% bar(0:.01:1, p_distribution)
% title('Distribution of P Values')
% xlabel('P value')
% ylabel('Frequency')
%--------------------------------------------------------------------------

%Find the ratio of the velocity component to the absolute value of the acceleration component. 
%Also find the ratio of the DFT as a function of direction.

unique_azimuth = unique(azimuth);
unique_elevation = unique(elevation);
b_to_a = vel_comp_b./(acc_comp_a);

k=1;
for i = -20:.1:20
    bbb(k) = length(find(b_to_a >= i & b_to_a < i+.1));
    k = k+1;
end


for j=1:length(unique_elevation)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) );
        b_to_a_ratio{i,j} = b_to_a(select);
        tuned_dir{i,j} = DFT_ratio(select);
    end
end
%--------------------------------------------------------------------------

%Find the mean for each direction

for j=1:length(unique_elevation)
    for i=1: length(unique_azimuth)
        if (j == 5 | j == 1) & (i ~=1)
            do_nothing = 1;
        else 
            %take the median to ignore outliers. outliers occur when a is
            %almost equal to zero.
            median_ba_ratio{i,j} = median(b_to_a_ratio{i,j});
            mean_tuned_dir{i,j} = mean(tuned_dir{i,j});  
        end   
    end
end

%first we throw out directions which have VAF's of less than 60% in the
%combined model.

directions = find(VAF_com >= 0.6);
good_dir = name(directions);

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
good_elevation = elevation(directions);
good_tau_com = tau_com(directions);
good_sig_com = sig_com(directions);


length(good_cells)
%The max direction will be the direction which has the highest DFT ratio.
for i = 1:length(good_cells)
    max_ratio = 0;
    for j = 1:length(good_dir)
        if (length(good_cells{i}) == length(good_dir{j})) & (sum(good_cells{i} == good_dir{j}) == length(good_cells{i})) & (max_ratio < new_DFT_ratio(j))
            max_ratio = new_DFT_ratio(j);
            max_az = good_azimuth(j);
            max_el = good_elevation(j);
            max_index(i) = j;
        end
    end
    %max_dir contains the max [Azimuth, Elevation] for each of the good
    %cells. Its index corresponds to that of good_cells such that
    %the max direction of good_cells{1} is max_dir{1}.
    
    max_dir{i} = [max_az, max_el];
    max_azimuth(i) = max_az;
    max_elevation(i) = max_el;
end
tau_max_com = good_tau_com(max_index);
%find angle between preferred direction and other directions
Azims = [0 45 90 135 180 225 270 315]*pi/180;
Elevs = [-90 -45 0 45 90]*pi/180;

%b/a ratio of each of the good directions
ba_vals = b_to_a(directions);

angle =[];
diff_angle{length(good_cells)} = [];
ba_rat{length(good_cells)} = [];

direc = zeros(0, length(good_dir));
count =1;

for i = 1:length(good_cells)
    for j = 1:length(good_dir)
        if (length(good_cells{i}) == length(good_dir{j})) & (sum(good_cells{i} == good_dir{j}) == length(good_cells{i}))
            temp_max_az = Azims(max_azimuth(i));
            temp_max_el = Elevs(max_elevation(i));
            temp_azimuth = Azims(good_azimuth(j));
            temp_elevation = Elevs(good_elevation(j));
            direc(count) = j;
            count = count+1;
            % calculate difference angles using the formula:
            % cos(x) = sin(Ea)sin(Eb) + cos(Eb)sin(Ab)cos(Ea)sin(Aa) + cos(Eb)cos(Ab)cos(Ea)cos(Aa)
            % where Ea and Eb is the elevation of vectors a and b, and Aa and Ab =
%             azimuth of a and b and x is the difference angle
            angle = round(acos(sin(temp_max_el)*sin(temp_elevation)+cos(temp_elevation)*sin(temp_azimuth)*cos(temp_max_el)*sin(temp_max_az)+cos(temp_elevation)*cos(temp_azimuth)*cos(temp_max_el)*cos(temp_max_az))*180/pi);
            
            diff_angle{i} = [diff_angle{i} angle];
            ba_rat{i} = [ba_rat{i} ba_vals(j)];            
        end
    end
end



for i = 1:length(diff_angle)
    unique_angle{i} = unique(diff_angle{i});
    
    for j = 1:length(unique_angle{i})
        find_index = find(diff_angle{i} == unique_angle{i}(j));
        good_ba_ratio{i}(j) = median(ba_rat{i}(find_index));
    end
end

%remove elements of unique anlge that that have only 2 data points since we
%cannot fit these.

% for i = 1:2
%     un5{i} = unique_angle{i};
%     ba5{i} = good_ba_ratio{i};
%     d_ang5{i} = diff_angle{i};
%     rat_ba5{i} = ba_rat{i};
% end
% for i = 4:84
%     un{i-3} = unique_angle{i};
%     ba{i-3} = good_ba_ratio{i};
%     d_ang{i-3} = diff_angle{i};
%     rat_ba{i-3} = ba_rat{i};
% end
% for i = 86:100
%     un2{i-85} = unique_angle{i};
%     ba2{i-85} = good_ba_ratio{i};
%     d_ang2{i-85} = diff_angle{i};
%     rat_ba2{i-85} = ba_rat{i};
% end
% for i = 102:107
%     un3{i-101} = unique_angle{i};
%     ba3{i-101} = good_ba_ratio{i};
%     d_ang3{i-101} = diff_angle{i};
%     rat_ba3{i-101} = ba_rat{i};
% end
% un4{1} = unique_angle{109};
% ba4{1} = good_ba_ratio{109};
% d_ang4{1} = diff_angle{109};
% rat_ba4{1} = ba_rat{109};
%     
% for i = 111:112
%     un6{i-110} = unique_angle{i};
%     ba6{i-110} = good_ba_ratio{i};
%     d_ang6{i-110} = diff_angle{i};
%     rat_ba6{i-110} = ba_rat{i};
% end
% 
% new_angle = [un5 un un2 un3 un4 un6];
% good_ba = [ba5 ba ba2 ba3 ba4 ba6];
% difference = [d_ang5 d_ang d_ang2 d_ang3 d_ang4 d_ang6];
% ba_points = [rat_ba5 rat_ba rat_ba2 rat_ba3 rat_ba4 rat_ba6];
% 
% difference{108} = [1 2 3];
% ba_points{108} = [3 2 1];

% for i = 1:6:103
% 
% %plot the change in b/a ratio as a function of the departure from the max direction with the max direction being 0   
%     
% %     figure
% %     plot(new_angle{i}, good_ba{i}, 'b', 'linewidth', 2)
% %     hold on
% %     plot(new_angle{i+1}, good_ba{i+1}, 'c', 'linewidth', 2)
% %     plot(new_angle{i+2}, good_ba{i+2}, 'g', 'linewidth', 2)
% %     plot(new_angle{i+3}, good_ba{i+3}, 'k', 'linewidth', 2)
% %     plot(new_angle{i+4}, good_ba{i+4}, 'm', 'linewidth', 2)
% %     strr = 'b/a ratio plot for cells:';
% %     to = '-';
% %     title([strr num2str(i) to num2str(i+4)])
% %     xlabel('Angle departing from pref. direction')
% %     ylabel('b/a ratio')
% %     orient landscape
% %     print -dwinc
% %     close
%         
% %     [b{i}, stats{i}] = robustfit(new_angle{i}', good_ba{i});
% %     y{1} = b{i}(2)*new_angle{i} + b{i}(1);
% %     [b{i+1}, stats{i+1}] = robustfit(new_angle{i+1}', good_ba{i+1});
% %     y{2} = b{i+1}(2)*new_angle{i+1} + b{i+1}(1);
% %     [b{i+2}, stats{i+2}] = robustfit(new_angle{i+2}', good_ba{i+2});
% %     y{3} = b{i+2}(2)*new_angle{i+2} + b{i+2}(1);
% %     [b{i+3}, stats{i+3}] = robustfit(new_angle{i+3}', good_ba{i+3});
% %     y{4} = b{i+3}(2)*new_angle{i+3} + b{i+3}(1);
% %     [b{i+4}, stats{i+4}] = robustfit(new_angle{i+4}', good_ba{i+4});
% %     y{5} = b{i+4}(2)*new_angle{i+4} + b{i+4}(1);
%     
% 
% %     
% %     figure
% %     plot(new_angle{i}, y{1}, 'b', 'linewidth', 2)
% %     hold on
% %     plot(new_angle{i+1}, y{2}, 'c', 'linewidth', 2)
% %     plot(new_angle{i+2}, y{3}, 'g', 'linewidth', 2)
% %     plot(new_angle{i+3}, y{4}, 'k', 'linewidth', 2)
% %     plot(new_angle{i+4}, y{5}, 'm', 'linewidth', 2)
% %     strr = 'Robust linear regression of b/a ratio for cells:';
% %     to = '-';
% %     title([strr num2str(i) to num2str(i+4)])
% %     xlabel('Angle departing from pref. direction')
% %     ylabel('b/a ratio')
% %     orient landscape
% %     print -dwinc
% %     close
% 
%     [b{i}, stats{i}] = robustfit(difference{i}', ba_points{i});
%     y{1} = b{i}(2)*difference{i} + b{i}(1);
%     [b{i+1}, stats{i+1}] = robustfit(difference{i+1}', ba_points{i+1});
%     y{2} = b{i+1}(2)*difference{i+1} + b{i+1}(1);
%     [b{i+2}, stats{i+2}] = robustfit(difference{i+2}', ba_points{i+2});
%     y{3} = b{i+2}(2)*difference{i+2} + b{i+2}(1);
%     [b{i+3}, stats{i+3}] = robustfit(difference{i+3}', ba_points{i+3});
%     y{4} = b{i+3}(2)*difference{i+3} + b{i+3}(1);
%     [b{i+4}, stats{i+4}] = robustfit(difference{i+4}', ba_points{i+4});
%     y{5} = b{i+4}(2)*difference{i+4} + b{i+4}(1);
%     [b{i+5}, stats{i+5}] = robustfit(difference{i+5}', ba_points{i+5});
%     y{6} = b{i+5}(2)*difference{i+5} + b{i+5}(1);
%     
% %     figure
% %     h(1) = plot(difference{i}, y{1}, 'b', 'linewidth', 2);
% %     hold on
% %     plot(difference{i}, ba_points{i}, 'ob', 'MarkerFaceColor', 'b')
% %     
% %     h(2) = plot(difference{i+1}, y{2}, 'c', 'linewidth', 2);
% %     plot(difference{i+1}, ba_points{i+1}, '^c', 'MarkerFaceColor', 'c')
% %     
% %     h(3) = plot(difference{i+2}, y{3}, 'm', 'linewidth', 2);
% %     plot(difference{i+2}, ba_points{i+2}, 'sm', 'MarkerFaceColor', 'm')
% %     
% %     h(4) = plot(difference{i+3}, y{4}, 'k', 'linewidth', 2);
% %     plot(difference{i+3}, ba_points{i+3}, 'dk', 'MarkerFaceColor', 'k')
% %     
% %     h(5) = plot(difference{i+4}, y{5}, 'r', 'linewidth', 2);
% %     plot(difference{i+4}, ba_points{i+4}, 'pr', 'MarkerFaceColor', 'r')
% %     
% %     h(6) = plot(difference{i+5}, y{6}, 'g', 'linewidth', 2);
% %     plot(difference{i+5}, ba_points{i+5}, 'hg', 'MarkerFaceColor', 'g')
% %        
% %     strr = 'Robust linear regression of b/a ratio for cells:';
% %     to = '-';
% %     title([strr num2str(i) to num2str(i+5)])
% %     xlabel('Angle departing from pref. direction')
% %     ylabel('b/a ratio')
% %     
% %     legend(h, num2str(stats{i}.p(2)), num2str(stats{i+1}.p(2)), num2str(stats{i+2}.p(2)), num2str(stats{i+3}.p(2)), num2str(stats{i+4}.p(2)), num2str(stats{i+5}.p(2)), 0)
%     
% %     orient landscape
% %     print -dwinc
% %     close
% 
%     
% end
% %--------------------------------------------------------------------------
% 
% for i = 1:107
%     params{i} = b{i};
%     statistics{i} = stats{i};
%     slopes(i) = params{i}(2);
%     regression_p_vals(i) = statistics{i}.p(2);
% end
% 
% %Plot a distribution of the slopes of all the cells.
% % tried doing this in a loop and stupid matlab keeps missing the zeros 
% bin_slope(1) = length(find(slopes > -.135 & slopes <= -.125));
% bin_slope(2) = length(find(slopes > -.125 & slopes <= -.115));
% bin_slope(3) = length(find(slopes > -.115 & slopes <= -.105));
% bin_slope(4) = length(find(slopes > -.105 & slopes <= -.095));
% bin_slope(5) = length(find(slopes > -.095 & slopes <= -.085));
% bin_slope(6) = length(find(slopes > -.085 & slopes <= -.075));
% bin_slope(7) = length(find(slopes > -.075 & slopes <= -.065));
% bin_slope(8) = length(find(slopes > -.065 & slopes <= -.055));
% bin_slope(9) = length(find(slopes > -.055 & slopes <= -.045));
% bin_slope(10) = length(find(slopes > -.045 & slopes <= -.035));
% bin_slope(11) = length(find(slopes > -.035 & slopes <= -.025));
% bin_slope(12) = length(find(slopes > -.025 & slopes <= -.015));
% bin_slope(13) = length(find(slopes > -.015 & slopes <= -.005));
% bin_slope(14) = length(find(slopes > -.005 & slopes <= .005));
% bin_slope(15) = length(find(slopes > .005 & slopes <= .015));
% bin_slope(16) = length(find(slopes > .015 & slopes <= .025));
% bin_slope(17) = length(find(slopes > .025 & slopes <= .035));
% bin_slope(18) = length(find(slopes > .035 & slopes <= .045));
% bin_slope(19) = length(find(slopes > .045 & slopes <= .055));
% bin_slope(20) = length(find(slopes > .055 & slopes <= .065));
% bin_slope(21) = length(find(slopes > .065 & slopes <= .075));
% bin_slope(22) = length(find(slopes > .075 & slopes <= .085));

% % figure
% % bar(-.13:.01:.08, bin_slope)
% % set(gca, 'xtick', [-.12 -.10 -.08 -.06 -.04 -.02 0 .02 .04 .06 .08])
% % xlabel('slopes')
% % ylabel('frequency')
% % title('Frequency Distribution of Slopes')
% %--------------------------------------------------------------------------
% 
% k = 1;
% for i = 0:.01:1
%     bin_p(k) = length(find((regression_p_vals <= i+.01) & (regression_p_vals > i)));
%     k = k+1;
% end

     %Scatter plot
good_VAF_com = VAF_com(directions);
good_VAF_vel = VAF_vel(directions);
%good_com_dir are the VAFs for the combined model directions that
%passed both the VAF and DFT criteria
good_veloc_dir = good_VAF_vel(direc);
good_com_dir = good_VAF_com(direc);

good_p1 = p(directions);
good_p = good_p1(direc);
significant = find(good_p <= 0.05);

% figure
% dum_vel2 = good_com_dir;
% plot(dum_vel2, good_com_dir, 'y', 'linewidth', 2)
% 
% hold on
% plot(good_veloc_dir, good_com_dir, 'bx') 
% xlabel('VAF velocity')
% ylabel('VAF combined')
% plot(good_veloc_dir(significant), good_com_dir(significant), 'rx')
% title('Scatter plor for VAFs for all directions with VAF > 60')
% orient landscape
% print -dwinc
% close


signif2 = find(p <= 0.05);

% figure
% dum_vel = VAF_com;
% plot(dum_vel, VAF_com, 'y', 'linewidth', 2)
% 
% hold on
% plot(VAF_vel, VAF_com, 'bx')
% xlabel('VAF velocity')
% ylabel('VAF combined')
% plot(VAF_vel(signif2), VAF_com(signif2), 'rx')
% title('Scatter plot of VAFs for all directions')
% orient landscape
% print -dwinc
% close


%--------------------------------------------------------------------------
%plot of mean b/a ratio vs the std of the b/a ratio
for i = 1:length(good_ba_ratio)
    median_ba(i) = median(good_ba_ratio{i});
    std_ba(i) = std(good_ba_ratio{i});
end

% figure
% plot(median_ba, std_ba, 'x')
% xlabel('Median b/a ratio')
% ylabel('STD of b/a ratio')
% title('STD of b/a ratio vs median of b/a ratio')
% axis([-1 (max(median_ba) + 1) -1 (max(std_ba) + 1)])

%--------------------------------------------------------------------------

%Find the matches of the good cells with Yong's cells; the cells which
%passed the VAF, 5 directions and DFT criteria

max_match = zeros(1,length(yong_cell_name));
for i = 1:length(yong_cell_name)
    for j = 1:length(good_cells)
        if (length(yong_cell_name{i}) == length(good_cells{j})) & (sum(yong_cell_name{i} == good_cells{j}) == length(yong_cell_name{i}))
            max_match(i) = 1;
        end
    end
end
max_match = find(max_match == 1);

yong_max_matched_cells = yong_cell_name(max_match);

%good cells matched with the pref dirs from the DFT method and yongs
%method. The pref direction is different from the max direction. The max
%direction is the direction that was recorded from which had a maximal
%response. The preferred direction is the direction calculated from the
%vector sum method. 
DFT_pref_matched = pref_cell_name(max_match);
DFT_pref_az = pref_az(max_match);
DFT_pref_el = pref_el(max_match);

yong_pref_matched = yong_pref_cell_name(max_match);
yong_pref_az = yong_pref_az(max_match);
yong_pref_el = yong_pref_el(max_match);

%max directions from yongs method
yong_max_el = yong_max_el(max_match);
yong_max_az = yong_max_az(max_match);
%--------------------------------------------------------------------------

%scatter plot that Dora wanted with just the max directions (for this to
%work, remember to select all cells that have AT LEAST 1 direction fitted.
%Not AT LEAST 5 like we have been doing for everything else). The figure
%from this is saved as scatterplot(max directions).fig so you dont have to
%run this part at all.

max_veloc_dir = good_VAF_vel(max_index);
max_com_dir = good_VAF_com(max_index);
max_p = good_p1(max_index);
signifi = find(max_p <= 0.05);
not_signifi = find(max_p > 0.05);
%find good cells which are significant with p<.05
significant_cells = find(max_p<=0.05);
% figure
% dum_vel2 = max_com_dir;
% plot(dum_vel2, max_com_dir, 'y', 'linewidth', 2)
% 
% hold on
% plot(max_veloc_dir, max_com_dir, 'bx') 
% xlabel('VAF velocity')
% ylabel('VAF combined')
% plot(max_veloc_dir(signifi), max_com_dir(signifi), 'rx')      
% title('Scatter plot max directions')
% orient landscape
% print -dwinc
% close
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
%Null = significant cells with free tau <= 1

%find the tau's for the max_directions
free_tau = tau_velocity(directions);
free_tau = free_tau(max_index);

k1=1;k2=1;k3=1;k4=1;k5=1;k6=1;
for i = 1:length(max_index)
    if (free_tau(i) > 1 & free_tau(i) < 2) & max_p(i) <= 0.05
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
    elseif free_tau(i) <= 1
        N(k5) = i;
        k5 = k5+1;
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%plot b/a ratio vs tau for max directions, separate by population
badum = b_to_a(directions);
b_a_max = badum(max_index);
% 
% figure
% plot(b_a_max, b_a_max, 'y', 'linewidth', 2)
% hold on
% plot(tau_max_com(A1), b_a_max(A1),'ro', 'markerfacecolor', 'r')
% plot(tau_max_com(A2), b_a_max(A2), 'go', 'markerfacecolor', 'g')
% plot(tau_max_com(A3), b_a_max(A3), 'ko', 'markerfacecolor', 'k')
% plot(tau_max_com(V), b_a_max(V), 'bo', 'markerfacecolor', 'b')
% %comment out this line for vetibular condition since there are not cells
% %with tau < 1
% plot(tau_max_com(N), b_a_max(N), 'mo', 'markerfacecolor', 'm')
% axis([1.5 3 -20 20])
% xlabel('Tau')
% ylabel('b/a ratio')
% title('b/a ratio vs. tau(combined) for max directions')
% legend('Unity' , 'A1', 'A2', 'A3', 'V', 'N')



%plot of figure from yong's paper of pref azimuth vs pref elevation

figure
plot(DFT_pref_az, DFT_pref_el, 'x')
hold on
plot(DFT_pref_az(signifi), DFT_pref_el(signifi), 'ro')
title('Plot of Preferred Azimuth vs Preferred Elevation with points with red circles indicating fits with p < 0.05')
xlabel('Preferred Azimuth')
ylabel('Preferred Elevation')
% orient landscape
% print -dwinc
% close
%--------------------------------------------------------------------------

%Distribution of preferred azimuths and elevations

k = 1;
for i = 0:30:360
    binned_total_azim(k) =  length(find(DFT_pref_az(not_signifi) >=i & DFT_pref_az(not_signifi) < i+30));
    binned_signif_azim(k) = length(find(DFT_pref_az(signifi) >=i & DFT_pref_az(signifi) < i+30));
k = k+1;
end

k = 1;
for i = -90:20:90
    binned_total_elev(k) = length(find(DFT_pref_el(not_signifi) >= i & DFT_pref_el(not_signifi) < i+20));
    binned_signif_elev(k) = length(find(DFT_pref_el(signifi) >=i & DFT_pref_el(signifi) < i+20));
    k = k+1;
end

figure
bar(7.5:30:367.5, binned_total_azim,0.5)
hold on 
bar(15+7.5:30:375+7.5, binned_signif_azim, 0.5,'r')
title('Distribution of the preferred azimuths (With significant cells shown in red)')
xlabel('Preferred Azimuth')
ylabel('Frequency')
set(gca,'xtick', [0 45 90 135 180 225 270 315 360])
legend('Non-significant', 'Significant')
% orient landscape
% print -dwinc
% close

% figure
% bar(-85:20:95, binned_total_elev,0.5)
% hold on 
% bar(-75:20:105, binned_signif_elev, 0.5,'r')
% title('Distribution of the preferred elevations (With significant cells shown in red)')
% xlabel('Preferred Elevation')
% ylabel('Frequency')
% set(gca, 'xtick', [-90 -45 0 45 90])
% legend('Non-significant', 'Significant')
% orient landscape
% print -dwinc
% close
%--------------------------------------------------------------------------
%plot the mean b/a ratio against preferred azimuth and elevation from the
%combined method
for i = 1:length(ba_rat)
    temp_ba_rat = ba_rat{i}(ba_rat{i} < 25);
mean_ba(i) = mean(temp_ba_rat);
end
max_mean_ba = mean_ba(significant_cells);
%wrap

DFT_pref_az(DFT_pref_az >= 270) = DFT_pref_az(DFT_pref_az >= 270) - 360;        
 
% figure
% plot(DFT_pref_az, mean_ba, 'bo')
% hold on
% plot(DFT_pref_az(significant_cells), max_mean_ba, 'bo', 'markerfacecolor', 'r')
% xlabel('Preferred Azimuth')
% ylabel('Mean b/a ratio')
% title('Mean b/a ratio vs. preferred azimuth')
% axis([-90 270 -10 10])
% set(gca,'xtick', [-90 0 45 90 135 180 225 270])
% legend('Non-significant cells', 'Significant cells',2)
% orient landscape
% print -dwinc
% close

% figure
% plot(DFT_pref_el, mean_ba, 'bo')
% hold on
% plot(DFT_pref_el(significant_cells), max_mean_ba, 'bo', 'markerfacecolor', 'r')
% xlabel('Preferred Elevation')
% ylabel('Mean b/a ratio')
% title('Mean b/a ratio vs. preferred elevation')
% axis([-90 90 -10 10])
% set(gca, 'xtick', [-90 -45 0 45 90])
% legend('Non-significant cells', 'Significant cells',2)
% orient landscape
% print -dwinc
% close
DFT_pref_az = pref_az(max_match);
%--------------------------------------------------------------------------

%plot distribution of the max directions 
az_dist = zeros(1,length(unique_azimuth));
el_dist = zeros(1,length(unique_elevation));
k=1;
for i = 0:45:360
    az_dist(k) = length(find(DFT_pref_az >= i & DFT_pref_az < i+45));
    k=k+1;
end

k=1;
for i = -90:45:90
    el_dist(k) = length(find(DFT_pref_el >= i & DFT_pref_el < i+45));
    k=k+1;
end

figure
bar(22.5:45:405, az_dist)
xlabel('Preferred Azimuth')
ylabel('Frequency')
title('Distribution of Azimuths')
set(gca,'xtick', [0 45 90 135 180 225 270 315 360])
% orient landscape
% print -dwinc
% close
% 
% figure
% bar(-67.5:45:135, el_dist)
% xlabel('Preferred Elevation')
% ylabel('Frequency')
% title('Distribution of Elevations')
% set(gca, 'xtick', [-90 -45 0 45 90])
% orient landscape
% print -dwinc
% close
%-------------------------------------------------------------------------------------            
%make a scatter plot of sigma and tau for the max directions only

tau_max_dir = good_tau_com(max_index);
sig_max_dir = good_sig_com(max_index);
dum_s = tau_max_dir;


figure
plot(tau_max_dir(not_signifi), sig_max_dir(not_signifi), 'x')
hold on
plot(tau_max_dir(signifi), sig_max_dir(signifi), 'rx')
title('Sigma vs Tau for the max directions')
xlabel('tau')
ylabel('sigma')
axis([0 3 0 1])
legend('Non-significant', 'Significant')
orient landscape
print -dwinc
close
%--------------------------------------------------------------------------

%distribution of tau and sigma

k=1;
for i = 2:.05:2.5
    signif_tau_dist(k) = length(find((tau_max_dir(signifi) >= i) & (tau_max_dir(signifi) < i + 0.05)));
    notsignif_tau_dist(k) = length(find((tau_max_dir(not_signifi) >= i) & (tau_max_dir(not_signifi) < i + 0.05)));
    k = k+1;
end
    
k = 1;
for i = 0:.05:max(sig_max_dir)
    signif_sig_dist(k) = length(find((sig_max_dir(signifi) >= i) & (sig_max_dir(signifi) < i + 0.05)));
    notsignif_sig_dist(k) = length(find((sig_max_dir(not_signifi) >= i) & (sig_max_dir(not_signifi) < i + 0.05)));
    k = k+1;
end

good_tau_vel = tau_velocity(directions);
tau_max_vel = good_tau_vel(max_index);
ss = find(p<= 0.05);
nn = find(p> 0.05);
max_sig = find(max_p <= 0.05);
max_not = find(max_p > 0.05);
%for cells fitted with only velocity model and tau unrestrained
k=1;
for i = 0:.05:3
    signif_tau_vel(k) = length(find((tau_velocity(ss) >= i) & (tau_velocity(ss) < i + 0.05)));
    notsignif_tau_vel(k) = length(find((tau_velocity(nn) >= i) & (tau_velocity(nn) < i + 0.05)));
    sig_max_tau(k) = length(find((tau_max_vel(max_sig) >= i) & (tau_max_vel(max_sig) < i + 0.05)));
    ns_max_tau(k) = length(find((tau_max_vel(max_not) >= i) & (tau_max_vel(max_not) < i + 0.05)));
    k = k+1;
end

% figure
% bar(2:0.05:2.5, notsignif_tau_dist,.5)
% hold on
% bar(2.025:.05:2.525, signif_tau_dist,.5,'r')
% xlabel('Tau')
% ylabel('Frequency')
% title('Distribution of Tau for max directions')
% % orient landscape
% % print -dwinc
% % close
% 
figure
bar(-0.025:.05:max(sig_max_dir)-.025, notsignif_sig_dist,0.5)
hold on
bar(0:.05:max(sig_max_dir), signif_sig_dist,0.5,'r')
xlabel('Sigma')
ylabel('Frequency')
title('Distribution of Sigma for max directions')
% orient landscape
% print -dwinc
% close

figure
bar(0:.05:3, notsignif_tau_vel,.5)
hold on
bar(.025:.05:3.025, signif_tau_vel,.5,'r')
xlabel('Tau')
ylabel('Frequency')
title('Distribution of Tau for all directions for just velocity')

figure
bar(0:.05:3, ns_max_tau,.5)
hold on
bar(.025:.05:3.025, sig_max_tau,.5,'r')
xlabel('Tau')
ylabel('Frequency')
title('Distribution of Tau for max directions for just velocity')

%--------------------------------------------------------------------------
% plot scatter plot of DFT preferred direction vs Yongs preferred direction

% % wrap data
% for i = 1:length(DFT_pref_az)
%     if DFT_pref_az(i) - yong_pref_az(i) > 180 & DFT_pref_az(i) > yong_pref_az(i)
%         DFT_pref_az(i) = DFT_pref_az(i) - 360;
%     end
%     if yong_pref_az(i) - DFT_pref_az(i) > 180 & yong_pref_az(i) > DFT_pref_az(i)
%         yong_pref_az(i) = yong_pref_az(i) - 360;
%     end
% end
    
% figure
% plot(yong_pref_az, yong_pref_az, 'y', 'linewidth', 2)
% hold on
% plot(yong_pref_az(not_signifi), DFT_pref_az(not_signifi), 'ro')
% plot(yong_pref_az(signifi), DFT_pref_az(signifi), 'co','markerfacecolor' , 'c')
% title('Scatter plot of Preferred Azimuths from DFT and Mean Firing rate methods')
% xlabel('Pref Azimuth (Mean Firing rate method)')
% ylabel('Pref Azimuth (DFT method)')
% axis([-135 360 -135 360])
% set(gca, 'xtick', -180:45:360)
% set(gca, 'ytick', -180:45:360)
% legend('Unity', 'Non-significant', 'Significant',2)
% orient landscape
% print -dwinc
% close

% figure
% plot(yong_pref_el, yong_pref_el, 'y', 'linewidth', 2)
% hold on
% plot(yong_pref_el(not_signifi), DFT_pref_el(not_signifi), 'ro')
% plot(yong_pref_el(signifi), DFT_pref_el(signifi), 'co', 'markerfacecolor' , 'c')
% title('Scatter plot of Preferred Elevations from DFT and Mean Firing rate methods')
% xlabel('Pref Elevation (Mean Firing rate method)')
% ylabel('Pref Elevation (DFT method)')
% axis([-90 90 -90 90])
% set(gca, 'xtick', -90:45:90)
% set(gca, 'ytick', -90:45:90)
% legend('Unity','Non-significant', 'Significant',2)
% orient landscape
% print -dwinc
% close
%--------------------------------------------------------------------------
%find the angle between mean firing rate pref dir and dft pref dir

for i = 1:length(max_match)
    % calculate difference angles using the formula:
    % cos(x) = sin(Ea)sin(Eb) + cos(Eb)sin(Ab)cos(Ea)sin(Aa) + cos(Eb)cos(Ab)cos(Ea)cos(Aa)
    % where Ea and Eb is the elevation of vectors a and b, and Aa and Ab =
    % azimuth of a and b and x is the difference angle
    Aa = yong_pref_az(i)*(pi/180);
    Ab = DFT_pref_az(i)*(pi/180);
    Ea = yong_pref_el(i)*(pi/180);
    Eb = DFT_pref_el(i)*(pi/180);
    direction_angle(i) = (acos(sin(Ea)*sin(Eb) + cos(Eb)*sin(Ab)*cos(Ea)*sin(Aa) + cos(Eb)*cos(Ab)*cos(Ea)*cos(Aa))*180/pi);
end
sig_direction_angle = direction_angle(signifi);
not_sig_direction_angle = direction_angle(not_signifi);
% plot a distribution of the delta preferred direction
k = 1;
for i = 0:15:180
    binned_dir_angle(k) = length(find(not_sig_direction_angle >= i & not_sig_direction_angle < i+15));
    binned_sig_dir_angle(k) = length(find(sig_direction_angle >= i & sig_direction_angle < i+15));
    k = k+1;
end

% figure
% bar(0:15:180, binned_dir_angle,.5)
% hold on
% bar(7.5:15:187.5, binned_sig_dir_angle, 0.5, 'r')
% title('Distribution of the difference between the Pref Directions from the Mean Firing Rate and DFT methods')
% xlabel('Difference in angle between Pref Directions')
% ylabel('Frequency')
% set(gca, 'xtick', 0:15:180)
% legend('Non-significant', 'Significant')
% orient landscape
% print -dwinc
% close

%--------------------------------------------------------------------------

%plot p vals for max directions
max_p_vals = good_p1(max_index);
k = 1;
for i = 0:0.05:.95
    p_max_dist(k) = length(find(max_p_vals >= i & max_p_vals < i+.05));
    k = k+1;
end

% figure
% bar(0.025:.05:.975, p_max_dist)
% title('Distribution of P values for max directions')
% xlabel('P values')
% ylabel('Frequency')
% axis([0 1 0 max(p_max_dist)])
% orient landscape
% print -dwinc
% close

%-------------------------------------------------------------------------
%plot direction_angle as a function of b/a ratio

max_direction_ba_ratio = ba_vals(max_index);
sig_max_dir_ba = max_direction_ba_ratio(signifi);
not_sig_max_dir_ba = max_direction_ba_ratio(not_signifi);

% figure
% plot(direction_angle(not_signifi), not_sig_max_dir_ba, 'x')
% hold on
% plot(direction_angle(signifi), sig_max_dir_ba, 'rx')
% title('Plot of difference angle vs. b/a ratio')
% xlabel('Difference Angle')
% ylabel('b/a ratio')
% axis([0 180 -20 20])
% legend('Non significant', 'Significant')
% orient landscape
% print -dwinc
% close

%_-------------------------------------------------------------------------

%also make a plot of the median b/a ratio over small ranges to see if the b/a
%ratio actually changes as a function of diff angle

k =1;
for i = 0:15:165
    dir_range{k} = find(direction_angle >= i & direction_angle < i+15);
    k = k+1;
end

for i = 1:length(dir_range)
    if length(dir_range{i}) > 0 
        median_ba_val(i) = median(max_direction_ba_ratio(dir_range{i}));
    else median_ba_val(i) = 0;
    end
end

% figure
% bar(7.5:15:172.5, median_ba_val)
% xlabel('Difference Angle')
% ylabel('Median b/a ratio')
% title('Plot of median b/a ratio vs difference angle')
% set(gca, 'xtick', 0:15:180)
% orient landscape
% print -dwinc
% close

%-------------------------------------------------------------------------
%get HTI figures
load HTI_data
k=1;
%serparate by significance of HTI's
for i = 0:.05:.6
    binned_vis_sig(k) = length(find(signif_vis>= i & signif_vis < i+.05));
    binned_vis_bad(k) = length(find(bad_vis>=i & bad_vis <i+.05));
    binned_vest_sig(k) = length(find(signif_vest>=i & signif_vest <i+.05));
    binned_vest_bad(k) = length(find(bad_vest >=i & bad_vest <i+.05));
    k=k+1;
end

HTI_vis =HTI_vis(max_match);
HTI_vest = HTI_vest(max_match);
%seperate by significance of fit model.
k =1;
for i = 0:.05:.6
    visual_sig_fit(k) = length(find(HTI_vis(signifi)>= i & HTI_vis(signifi) < i+.05));
    visual_not_sig_fit(k) = length(find(HTI_vis(not_signifi) >=i & HTI_vis(not_signifi) <i+.05));
    vest_sig_fit(k) = length(find(HTI_vest(signifi) >=i & HTI_vest(signifi) <i+.05));
    vest_not_sig_fit(k) = length(find(HTI_vest(not_signifi) >=i & HTI_vest(not_signifi) <i+.05));
    k=k+1;
end

% figure
% bar(0:.05:.6, binned_vis_sig,.5,'r')
% hold on
% bar(0.025:.05:.625,binned_vis_bad, .5) 
% axis([0 .6 0 max([max(binned_vis_bad) max(binned_vis_sig)])])
% xlabel('HTI')
% ylabel('Number of Neurons')
% title('Visual Condition')
% legend('p<.05', 'p>.05')
% orient landscape
% print -dwinc
% close


% figure
% bar(0:.05:.6, binned_vest_sig, .5,'r')
% hold on
% bar(0.025:.05:.625, binned_vest_bad, .5)
% axis([0 .6 0 max([max(binned_vest_bad) max(binned_vest_sig)])])
% xlabel('HTI')
% ylabel('Number of Neurons')
% title('Vestibular Condition')
% legend('p<.05', 'p>.05')
% orient landscape
% print -dwinc
% close

% figure
% bar(0.0125:.05:.6125, visual_sig_fit,.5,'r')
% hold on
% bar(0.0375:.05:.6375,visual_not_sig_fit, .5) 
% axis([0 .6 0 max([max(visual_not_sig_fit) max(visual_sig_fit)])])
% xlabel('HTI distribution for max directions')
% ylabel('Number of Neurons')
% title('HTI distribution for max directions: Visual Condition')
% legend('p<.05', 'p>.05')
% orient landscape
% print -dwinc
% close


% figure
% bar(0.0125:.05:.6125, vest_sig_fit, .5,'r')
% hold on
% bar(0.0375:.05:.6375, vest_not_sig_fit, .5)
% axis([0 .6 0 max([max(vest_not_sig_fit) max(vest_sig_fit)])])
% xlabel('HTI')
% ylabel('Number of Neurons')
% title('HTI distribution for max directions: Vestibular Condition')
% legend('p<.05', 'p>.05')
% orient landscape
% print -dwinc
% close


%--------------------------------------------------------------------------
%to use this part of the code, first run analyzer.m with the first few
%lines uncommented and the rest commented, to generate ba_vestibular and
%good_vestibular cells. Remember to change the data that is being loaded
%aswell

% -------------------------------------------------------------
% ba_vestibular = ba_vals(max_index);
% good_vestibular_cells = good_cells;
% vestibular_azimuth = azimuth(directions);
% vestibular_elevation = elevation(directions);
% max_index_vestibular = max_index;
% max_match_vestibular = max_match;
% free_tau_vest = free_tau;
% max_p_vest = max_p;
% load HTI_data
% index_vis = zeros(1,255);
% index_vest = zeros(1,255);
% index_vest(find(p_vest <= 0.05)) = 1;
% index_vis(find(p_vis <=0.05)) = 1;
% index = index_vest+index_vis;
% index = find(index == 2);
% vest_az = pref_az(index);
% vest_el = pref_el(index);
% direction_angle_vest = direction_angle;
%------------------------------------------------------------- 




% % %we need to find the cells that have been fitted in both the vestibular and
% % %visual conditions in the max direction which is defined as the visual max
% % %direction.
% 
% vis_vest_match = zeros(1, length(good_cells));
% for i = 1:length(good_vestibular_cells)
%     for j = 1:length(good_cells)
%         if (length(good_vestibular_cells{i}) == length(good_cells{j})) & (sum(good_vestibular_cells{i} == good_cells{j}) == length(good_vestibular_cells{i}))
%             vis_vest_match(j) = 1;
%         end
%     end
% end
% 
% matched_cell_index = find(vis_vest_match == 1);
% 
% dummy_cells = good_cells(matched_cell_index);
% for i = 1:length(dummy_cells)
%     for j = 1:length(good_vestibular_cells)
%         if (length(dummy_cells{i}) == length(good_vestibular_cells{j})) & (sum(dummy_cells{i} == good_vestibular_cells{j}) == length(dummy_cells{i}))
%             vis_vest_match2(j) = 1;
%         end
%     end
% end
% matched_cell_index2 = find(vis_vest_match2 == 1);
% 
% free_tau = free_tau(matched_cell_index);
% max_p = max_p(matched_cell_index);
% k1=1;k2=1;k3=1;k4=1;k5=1;k6=1;
% for i = 1:length(matched_cell_index)
%     if (free_tau(i) > 1 & free_tau(i) < 2) & max_p(i) <= 0.05
%         A1_com(k1) = i;
%         k1=k1+1;
%     elseif (free_tau(i) >= 2 & free_tau(i) < 2.25) & max_p(i) <= 0.05
%         A2L_com(k2) = i;
%         k2=k2+1;
%     elseif (free_tau(i) >= 2.25 & free_tau(i) < 2.5) & max_p(i) <= 0.05
%         A2R_com(k6) = i;
%         k6=k6+1;
%     elseif (free_tau(i) >= 2.5 & free_tau(i) <= 3) & max_p(i) <= 0.05
%         A3_com(k3) = i;
%         k3=k3+1;
%     elseif free_tau(i) > 1 & max_p(i) > 0.05
%         V_com(k4) = i;
%         k4=k4+1;
%     elseif free_tau(i) <= 1
%         N_com(k5) = i;
%         k5 = k5+1;
%     end
% end
% 
% ba_visual = ba_vals(max_index);
% ba_visual = ba_visual(matched_cell_index);
% 
% ba_vestibular = ba_vestibular(matched_cell_index2);
% % 
% % 
% % % figure
% % % plot(ba_visual, ba_visual, 'y', 'linewidth', 2)
% % % hold on
% % % plot(ba_visual(A1_com), ba_vestibular(A1_com), 'ro', 'markerfacecolor', 'r')
% % % plot(ba_visual(A2L_com), ba_vestibular(A2L_com), 'go', 'markerfacecolor', 'g')
% % % plot(ba_visual(A2R_com), ba_vestibular(A2R_com), 'co', 'markerfacecolor', 'c')
% % % plot(ba_visual(A3_com), ba_vestibular(A3_com), 'ko', 'markerfacecolor', 'k')
% % % plot(ba_visual(V_com), ba_vestibular(V_com), 'bo', 'markerfacecolor', 'b')
% % % plot(ba_visual(N_com), ba_vestibular(N_com), 'mo', 'markerfacecolor', 'm')
% % % axis([-15 15 -15 15])
% % % title('Scatter plot for b/a ratios for visual and vestibular conditions for all max directions')
% % % xlabel('Visual b/a ratio')
% % % ylabel('Vestibular b/a ratio')
% % % legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V', 'N')
% % %--------------------------------------------------------------------------
% % %scatter plot of yong HTI vs DFT HTI (uses the indexes from above so to run
% % %this you have to first run vestibular then visual as above.)
% % 
% load HTI_data
% 
% order_index = zeros(1,length(name));
% for i = 1:length(name)
%     for j = 1:length(name_yong)
%         if (length(name_yong{j}) == length(name{i})) & (sum(name_yong{j} == name{i}) == length(name_yong{j}))
%             order_index(i) = 1;
%         end
%     end
% end
% order_index = find(order_index == 1);
% HTI_vis_yong = HTI_vis_yong(order_index);
% HTI_vest_yong = HTI_vest_yong(order_index);
% 
% HTI_vis_yong = HTI_vis_yong(max_match);
% HTI_vis_yong = HTI_vis_yong(matched_cell_index);
% HTI_vis = HTI_vis(max_match);
% HTI_vis = HTI_vis(matched_cell_index);
% 
% HTI_vest_yong = HTI_vest_yong(max_match_vestibular);
% HTI_vest_yong = HTI_vest_yong(matched_cell_index2);
% HTI_vest = HTI_vest(max_match_vestibular);
% HTI_vest = HTI_vest(matched_cell_index2);
% 
% figure
% plot(HTI_vis_yong, HTI_vis_yong, 'y', 'linewidth', 2);
% hold on
% plot(HTI_vis(A1_com), HTI_vis_yong(A1_com), 'ro', 'markerfacecolor', 'r')
% plot(HTI_vest(A2L_com), HTI_vest_yong(A2L_com), 'go', 'markerfacecolor', 'g')
% plot(HTI_vest(A2R_com), HTI_vest_yong(A2R_com), 'co', 'markerfacecolor', 'c')
% plot(HTI_vis(A3_com), HTI_vis_yong(A3_com), 'ko', 'markerfacecolor', 'k')
% plot(HTI_vis(V_com), HTI_vis_yong(V_com), 'bo', 'markerfacecolor', 'b')
% % plot(HTI_vis(N_com), HTI_vis_yong(N_com), 'mo', 'markerfacecolor', 'm')
% axis([0 0.8 0 0.8])
% xlabel('DFT method HTI')
% ylabel('Mean Firing Rate method HTI')
% title('Scatter plot of HTIs from mean firing rate and DFT methods for max directions- VISUAL CONDITION')
% legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V')
% 
% 
% free_tau_vest = free_tau_vest(matched_cell_index2);
% max_p_vest = max_p_vest(matched_cell_index2);
% k1=1;k2=1;k3=1;k4=1;k5=1;k6=1;
% for i = 1:length(matched_cell_index2)
%     if (free_tau_vest(i) > 1 & free_tau_vest(i) < 2) & max_p_vest(i) <= 0.05
%         A1_vest(k1) = i;
%         k1=k1+1;
%     elseif (free_tau_vest(i) >= 2 & free_tau_vest(i) < 2.25) & max_p_vest(i) <= 0.05
%         A2L_vest(k2) = i;
%         k2=k2+1;
%     elseif (free_tau_vest(i) >= 2.25 & free_tau_vest(i) < 2.5) & max_p_vest(i) <= 0.05
%         A2R_vest(k6) = i;
%         k6=k6+1;
%     elseif (free_tau_vest(i) >= 2.5 & free_tau_vest(i) <= 3) & max_p_vest(i) <= 0.05
%         A3_vest(k3) = i;
%         k3=k3+1;
%     elseif free_tau_vest(i) > 1 & max_p_vest(i) > 0.05
%         V_vest(k4) = i;
%         k4=k4+1;
%     elseif free_tau_vest(i) <= 1
%         N_vest(k5) = i;
%         k5 = k5+1;
%     end
% end
% 
% figure
% plot(HTI_vest_yong, HTI_vest_yong, 'y', 'linewidth', 2);
% hold on
% plot(HTI_vest(A1_vest), HTI_vest_yong(A1_vest), 'ro', 'markerfacecolor', 'r')
% plot(HTI_vest(A2L_vest), HTI_vest_yong(A2L_vest), 'go', 'markerfacecolor', 'g')
% plot(HTI_vest(A2R_vest), HTI_vest_yong(A2R_vest), 'co', 'markerfacecolor', 'c')
% plot(HTI_vest(A3_vest), HTI_vest_yong(A3_vest), 'ko', 'markerfacecolor', 'k')
% plot(HTI_vest(V_vest), HTI_vest_yong(V_vest), 'bo', 'markerfacecolor', 'b')
% axis([0 0.8 0 0.8])
% xlabel('DFT method HTI')
% ylabel('Mean Firing Rate method HTI')
% title('Scatter plot of HTIs from mean firing rate and DFT methods for max directions- VESTIBULAR CONDITION')
% legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V')
% 
% %--------------------------------------------------------------------------
% %distribution of diff angle from vis/vest conditions
% 
% % visual_az = pref_az(index);
% % visual_el = pref_el(index);
% % for i = 1:length(visual_az)
% %     Aa = visual_az(i)*(pi/180);
% %     Ab = vest_az(i)*(pi/180);
% %     Ea = visual_el(i)*(pi/180);
% %     Eb = vest_el(i)*(pi/180);
% %     dif_ang_vis_vest(i) = (acos(sin(Ea)*sin(Eb) + cos(Eb)*sin(Ab)*cos(Ea)*sin(Aa) + cos(Eb)*cos(Ab)*cos(Ea)*cos(Aa))*180/pi);
% % end
% % 
% % % plot a distribution of the delta preferred direction
% % k = 1;
% % for i = 0:15:180
% %     binned_angle(k) = length(find(dif_ang_vis_vest >= i & dif_ang_vis_vest < i+15));
% %     k = k+1;
% % end
% % figure
% % bar(0:15:180, binned_angle)
% % title('Distribution of the difference between the Pref Directions from the visual and vestibular conditions')
% % xlabel('Difference in angle between conditions')
% % ylabel('Number')
% % set(gca, 'xtick', 0:15:180)
% %--------------------------------------------------------------------------
% %scatter plot of diff angle from vis/vest conditions
% 
% direction_angle_vis = direction_angle(matched_cell_index);
% direction_angle_vest = direction_angle_vest(matched_cell_index2);
% figure
% plot(direction_angle_vis, direction_angle_vis, 'y', 'linewidth', 2)
% hold on
% plot(direction_angle_vis(A1_com),direction_angle_vest(A1_com),'o', 'markerfacecolor', 'b')
% plot(direction_angle_vis(A2L_com),direction_angle_vest(A2L_com),'ro', 'markerfacecolor', 'r')
% plot(direction_angle_vis(A2R_com),direction_angle_vest(A2R_com),'go', 'markerfacecolor', 'g')
% plot(direction_angle_vis(A3_com),direction_angle_vest(A3_com),'k^', 'markerfacecolor', 'k')
% plot(direction_angle_vis(V_com),direction_angle_vest(V_com),'m^', 'markerfacecolor', 'm')
% legend('Unity', 'A1', 'A2L', 'A2R', 'A3', 'V')
% title('Scatter plot of diff angle btwn mean firing rate and DFT method- visual vs. vest')
% xlabel('Visual Difference Angle')
% ylabel('Vestibular Difference Angle')