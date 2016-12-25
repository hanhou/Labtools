clear all;

filepath = 'C:\MATLAB6p5\work\';
filename = 'PSTH_fit_ves.txt';

[FILE azimuth gaze p_val VAF_com VAF_vel DC_com b_com tau_com sigma_com a_com DC_vel K_vel tau_vel sigma_vel DC_unc K_unc tau_unc sigma_unc DFT_ratio] = textread([filepath filename], ... 
'%s %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', '\t', 'headerlines', 1);

file = cell2mat(FILE);
unique_file = unique(FILE);

index = 1;
for n = 1:length(unique_file)
    max_ratio = 0;
    while file(index,:) == cell2mat(unique_file(n))
        if DFT_ratio(index) > max_ratio
            prefdir_VAF_com(n) = VAF_com(index);
            prefdir_VAF_vel(n) = VAF_vel(index);
            prefdir_VAF_diff(n) = prefdir_VAF_com(index) - prefdir_VAF_vel(index);
            prefdir_freetau(n) = tau_unc(index);
            prefdir_b(n) = b_com(index);
            prefdir_a(n) = a_com(index);
            prefdir_ratio(n) = DFT_ratio(index);
            max_ratio = DFT_ratio(index);
        end
        index = index + 1;
        if index > length(file)
            break
        end
    end
end

prefdir_VAF_com = prefdir_VAF_com';
prefdir_VAF_vel = prefdir_VAF_vel';
prefdir_VAF_diff = prefdir_VAF_diff';
prefdir_freetau = prefdir_freetau';
prefdir_b = prefdir_b';
prefdir_a = prefdir_a';
prefdir_ratio = prefdir_ratio';


% Z_com = 0.5 * log( (1+sqrt(prefdir_VAF_com)) / (1-sqrt(prefdir_VAF_com)) );
% Z_vel = 0.5 * log( (1+sqrt(prefdir_VAF_vel)) / (1-sqrt(prefdir_VAF_vel)) );
% Z_com2 = 0.5 * ( log(1+sqrt(prefdir_VAF_com)) - log(1-sqrt(prefdir_VAF_com)) );