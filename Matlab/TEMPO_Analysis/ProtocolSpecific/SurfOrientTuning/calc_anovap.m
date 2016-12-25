function p_val = calc_anovap(plot_x, plot_y) 
spike(:,1) = plot_x';
spike(:,2) = plot_y';
sortedspikes = sortrows(spike, [1]);
p_val = anovan(sortedspikes(:, 2), {sortedspikes(:, 1)}, 'full', 3, {'Tilt'}, 'off');