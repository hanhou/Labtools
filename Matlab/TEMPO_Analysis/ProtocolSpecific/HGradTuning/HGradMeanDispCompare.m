function calc_mdisp_anovap(disp_select, plot_x, plot_y, unique_disp_ang) 
spike(:,1) = plot_x';
spike(:,2) = plot_y';
sortedspikes = sortrows(spike, [1]);
for temp_ang = 1:length(unique_disp_ang)
    ang_ind = find(sortedspikes(:,1) == unique_disp_ang(temp_ang));
    sortedspikes(ang_ind(1):ang_ind(length(ang_ind)), 1) = temp_ang;
end
p_val = anovan(sortedspikes(:, 2), {sortedspikes(:, 1)}, 'full', 3, {'Angles'}, 'off');
