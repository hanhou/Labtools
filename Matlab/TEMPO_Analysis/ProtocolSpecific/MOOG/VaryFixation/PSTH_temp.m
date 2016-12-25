clear all;

filepath = 'C:\MATLAB6p5\work\';
filename = 'PSTH_fit_vis.txt';

[FILE azimuth gaze p_val VAF_com VAF_vel DC_com b_com tau_com sigma_com a_com DC_vel K_vel tau_vel sigma_vel DC_unc K_unc tau_unc sigma_unc DFT_ratio] = textread([filepath filename], ... 
'%s %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', '\t', 'headerlines', 1);

file = cell2mat(FILE);
unique_file = unique(FILE);

index = 1;
for n = 1:length(unique_file)
    max_ratio = 0;
    while file(index,:) == cell2mat(unique_file(n)) & index <= length(FILE)
        if DFT_ratio(index) > max_ratio
            prefdir_VAF_com(n) = VAF_com(index);
            prefdir_VAF_diff(n) = VAF_com(index) - VAF_vel(index);
            prefdir_ratio(n) = DFT_ratio(index);
            max_ratio = DFT_ratio(index);
        end
        index = index + 1;
        if index > length(file)
            break
        end
    end
end

prefdir_ratio = prefdir_ratio';
prefdir_VAF_com = prefdir_VAF_com';
prefdir_VAF_diff = prefdir_VAF_diff';