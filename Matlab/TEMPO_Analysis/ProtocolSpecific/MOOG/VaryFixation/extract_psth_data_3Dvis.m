clear all;

filepath = 'Z:\Users\Anuk\';
filename = 'fit_data_chris.mat';

load([filepath filename]);

name{370} = 'm3c016r1.htb';
name{371} = 'm3c016r1.htb';
name{372} = 'm3c016r1.htb';
name{373} = 'm3c016r1.htb';
name{374} = 'm3c016r1.htb';
name{375} = 'm3c016r1.htb';
name{376} = 'm3c016r1.htb';
name{377} = 'm3c016r1.htb';

file = cell2mat(name);
unique_file = unique(name);

index = 1;
for n = 1:length(unique_file)
    max_ratio = 0;
    while file(index,:) == cell2mat(unique_file(n))
        if DFT_ratio(index) > max_ratio & elevation(index) == 3
            prefdir_VAF_com(n) = VAF_com(index);
            prefdir_VAF_vel(n) = VAF_vel(index);
            prefdir_VAF_diff(n) = prefdir_VAF_com(n) - prefdir_VAF_vel(n);
            prefdir_freetau(n) = tau_velocity(index);
            prefdir_b(n) = vel_comp_b(index);
            prefdir_a(n) = acc_comp_a(index);
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
