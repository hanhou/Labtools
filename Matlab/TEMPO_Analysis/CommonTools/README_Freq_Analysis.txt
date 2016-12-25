%%%%%%%%% EDITED BY CRF, 12/04/07 %%%%%%%%%%%%
Added null trial related vars, p vals for comparison of DFTR and null-trial DFTR, etc.,
and some fixes to the text, including new folder/data saving method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% EDITED BY CRF, 3/27/07 %%%%%%%%%%%%
Added a few variables to the output (marked with *** below), and cleaned up the text.
(NOTE: The output (especially the 'textread' code below) will not be compatible with 
pre-existing data files.  Backup your old data and start from scratch.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


To run this analysis, execute the "Run Frequency Analysis" protocol in the 3D Direction Tuning protocol in batch_gui or tempo_gui (if you have only one cell).  If it is the rotation protocol, select "Rotation Frequency Analysis".  Two data files will be generated in the path specified at the end of the "Frequency_Analysis.m" (or Rotation_Frequency_Analysis.m )file. 


The two data files will contain all the data that you need.  To read from the data files, copy and paste the text below (beginning with 'filedir = ...') into the command window and the data will be loaded into the workspace.  If your data files have column headers, you must remove them before running the commands below.  Also remove the 'word-wrap' formatting on this text file before copying:

***IMPORTANT: modify the 'filedir' string to a folder on your hard drive.***


filedir = 'C:\MATLAB6p5\work\Frequency_Analysis\versions\';
FILE1 = [filedir fileversion '_Frequency_data1.dat'];
FILE2 = [filedir fileversion '_Frequency_data2.dat'];

str1 = '%s %f %f %f %f %f %f %f %f %f %f %f %f %f';
[name, azimuth, elevation, dft_ratio, dft_p, dft_vel, dft_p_vel, dft_acc, dft_p_acc, phase, delay, p_dft_vs_null, p_vel_vs_null, p_acc_vs_null] = textread(FILE1, str1);
str2 = '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
[name2, P_anova, P_anova_full, P_anova_vel, P_anova_acc, DDI_mfr, DDI_full, DDI_vel, DDI_acc, p_ddi_mfr, p_ddi_full, p_ddi_vel, p_ddi_acc, corr_vel_mfr, p_vel_mfr, corr_acc_mfr, p_acc_mfr, corr_acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc] = textread(FILE2, str2);




----------------------------------------------------------------------------
The data in "Frequency_data1.dat" file contain data for each response direction for each cell. Thus each variable will contain 27 consecutive rows that correspond to the 26 motion directions tested, plus the null trials (designated 9999). Data are stored from directions (Azim, Elev) recorded in the following order:

   Azim		Elev
     0		90				
     0		45		
     0		0
     0		-45
     0		-90
    45		45
    45		0
    45		-45
    90		45		
    90		0
    90		-45
   135		45
   135		0
   135		-45
   180		45
   180		0
   180		-45
   225		45
   225		0
   225		-45
   270		45
   270		0
   270		-45
   315		45
   315		0
   315		-45
   9999		9999


14 VARIABLES IN Frequency_data1.dat:

name - cell name
*** azimuth - azimuth of direction or axis of rotation
*** elevation - elevation of direction or axis of rotation
dft_ratio - DFT ratio ('full')
dft_p - p value for DFT ratio (permutation test)
*** dft_vel - velocity component of DFTR 
dft_p_vel - p value for velocity component of DFTR (permutation test)
*** dft_acc - acceleration component of DFTR
dft_p_acc - p value for acceleration component of DFTR (permutation test)
phase - phase of response (largest component)
delay - delay of response (latency to peak)
p_dft_vs_null - p val for ANOVA of single-trial DFTR for each dir vs. DFTR of the null trials
p_vel_vs_null - same for velocity
p_acc_vs_null - same for acceleration

----------------------------------------------------------------------------
The data in "Frequency_data2.dat" file contain one data entry for each cell.

28 VARIABLES IN Frequency_data2.dat:

name2 - cell name
P_anova - p value from anova on mean firing rate (spatial tuning significance)
*** P_anova_full - p value from anova on 'full' DFTR
P_anova_vel - p value from anova on DFTR velocity component
P_anova_acc - p value from anova on DFTR acceleration component
*** DDI_mfr - DDI based on mean firing rate
*** DDI_full - DDI based on full DFTR
DDI_vel - DDI based on DFTR velocity component 
DDI_acc - DDI based on DFTR acceleration component 
*** p_ddi_mfr - p value for DDI_mfr (permutation test)
*** p_ddi_full - p value for DDI_full (permutation test) 
p_ddi_vel - p value for DDI_vel (permutation test)
p_ddi_acc - p value for DDI_acc (permutation test)
corr_vel_mfr - correltaion coefficient of 3D tuning curves for velocity vs mean firing rate
p_vel_mfr - p value for corr_vel_mfr
corr_acc_mfr - correltaion coefficient of 3D tuning curves for acceleration vs mean firing rate
p_acc_mfr - p value for corr_acc_mfr
corr_acc_vel - correltaion coefficient of 3D tuning curves for velocity vs acceleration
p_acc_vel - p value for corr_acc_vel 
mfr_pref_az - preferred azimuth based on mfr 
mfr_pref_el - preferred elevation based on mfr 
mfr_pref_amp - preferred amplitude based on mfr 
pref_az_vel - preferred azimuth based on DFTR velocity component
pref_el_vel - preferred elevation based on DFTR velocity component 
pref_amp_vel - preferred amplitude based on DFTR velocity component
pref_az_acc - preferred azimuth based on DFTR acceleration component 
pref_el_acc - preferred elevation based on DFTR acceleration component
pref_amp_acc - preferred amplitude based on DFTR acceleration component