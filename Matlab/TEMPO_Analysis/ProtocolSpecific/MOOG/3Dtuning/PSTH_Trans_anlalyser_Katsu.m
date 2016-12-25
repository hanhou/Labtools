datasize = dlmread('Peak_Trans_PSTH_output_Katsu.dat','',1,1);
[row col]=size(datasize) 

vest = dlmread('Peak_Trans_PSTH_output_Katsu.dat','',[1 1 row 40]);size(vest)
visu = dlmread('Peak_Trans_PSTH_output_Katsu.dat','',[1 41 row 80]);size(visu)
%  aa3 = dlmread('Peak_PSTH_output_Katsu.dat','',[1 81 1 120]);dim=size(aa3)
figure;bar(mean(vest))
figure;bar(mean(visu))