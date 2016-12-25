clear all;

filepath = 'Z:\Users\Anuk\';
filename = 'fit_data_vestibular_chris.mat';

load([filepath filename]);

name_chris = sort(name_chris);
file = cell2mat(name_chris);
unique_file = unique(name_chris);

index = 1;
for n = 1:length(unique_file)
    max_ratio = 0;
    while file(index,:) == cell2mat(unique_file(n))
        if DFT_ratio_chris(index) > max_ratio & elevation_chris(index) == 3
            prefdir_VAF_com(n) = VAF_com_chris(index);
            prefdir_VAF_vel(n) = VAF_vel_chris(index);
            prefdir_VAF_diff(n) = prefdir_VAF_com(n) - prefdir_VAF_vel(n);
            prefdir_freetau(n) = tau_velocity_chris(index);
            prefdir_b(n) = vel_comp_b_chris(index);
            prefdir_a(n) = acc_comp_a_chris(index);
            prefdir_ratio(n) = DFT_ratio_chris(index);
            prefdir_file(n,:) = file(index,:);
            max_ratio = DFT_ratio_chris(index);
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
