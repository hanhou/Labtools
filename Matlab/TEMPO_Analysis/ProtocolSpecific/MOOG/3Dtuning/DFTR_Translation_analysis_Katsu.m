% This m file plots various quantities obtained by running
% Frequency_analysis, Rotation_frequency_analysis or frequency_analysis_2d.
% written by Narayan Ganesan (2007)
% REMEMBER to change the textread filenames to your specific files.(some notes - ab)

clear all;
str1 =  ' %s%f%f%f%f%f%f'; %6 values for direction specific (26/cell for the rotation and direction tuning protocols)
[name  dft_ratio, dft_p, dft_p_vel, dft_p_acc, phase, delay] = textread('Translation_Frequency_26trajectories.dat', str1);
str1 =  ' %s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'; %22 values for cell specific
[name2 P_anova, P_anova_vel, P_anova_acc, DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, corr_vel_mfr, p_vel_mfr, corr_acc_mfr, p_acc_mfr, corr_acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc] = textread('Translation_Frequency_cell.dat', str1);

vel_comp = -dft_ratio.*cos(phase);
acc_comp = dft_ratio.*sin(phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   insert all 26 trajectories
cellno=length(dft_ratio)/26;%n = No. of cell, temporaly 46 cells (03/07/07)
for n=1:cellno;%n = No. of cell, temporaly 46 cells (03/07/07)
    for m=1:26;% m=trajectories
        dft_ratio26{n}(m)=dft_ratio((n-1)*26+m);
        dft_ratio26_vel{n}(m)=vel_comp((n-1)*26+m);
        dft_ratio26_acc{n}(m)=acc_comp((n-1)*26+m);
        dft_p26{n}(m)=dft_p((n-1)*26+m);
        dft_p26_vel{n}(m)=dft_p_vel((n-1)*26+m);
        dft_p26_acc{n}(m)=dft_p_acc((n-1)*26+m);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_sig_dft(n)=0;
count_sig_dft_vel(n)=0;
count_sig_dft_acc(n)=0;
for n=1:(length(dft_ratio)/26);%n = No. of cell, temporaly 46 cells (03/07/07)
    for m=1:26;% m=trajectories
        if dft_p26{n}(m)<0.05
            count_sig_dft(n)=count_sig_dft(n)+1;
        end
        if dft_p26_vel{n}(m)<0.05
            count_sig_dft_vel(n)=count_sig_dft_vel(n)+1;
        end
        if dft_p26_acc{n}(m)<0.05
            count_sig_dft_acc(n)=count_sig_dft_acc(n)+1;
        end
    end
%     count_sig_dft(n)
end
% P_anova=P_anova'
P_anova
count_sig_dft=count_sig_dft'
P_anova_vel
count_sig_dft_vel=count_sig_dft_vel'
P_anova_acc
count_sig_dft_acc=count_sig_dft_acc'
% now copy paste to origin is better!
%  figure(2);
% % semilogy(count_sig_dft,P_anova,'o');hold on;
% plot(count_sig_dft,P_anova,'o');hold on;
% startpoint=[0 26];endpoint=[0.05 0.05]; line(startpoint,endpoint);hold on;
% xlim([0,26]);ylim([0, 1]);
% 
% figure(3);
% % semilogy(count_sig_dft,P_anova,'o');hold on;
% plot(count_sig_dft_vel,P_anova,'o');hold on;
% startpoint=[0 26];endpoint=[0.05 0.05]; line(startpoint,endpoint);hold on;
% xlim([0,26]);ylim([0, 1]);
% 
% figure(4);
% % semilogy(count_sig_dft,P_anova,'o');hold on;
% plot(count_sig_dft_acc,P_anova,'o');hold on;
% startpoint=[0 26];endpoint=[0.05 0.05]; line(startpoint,endpoint);hold on;
% xlim([0,26]);ylim([0, 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% convert i, j
%%% change 26 trajectories to azix8 and elex5
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];
for n=1:cellno
% n=1;
for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                DFT_mfr{n}(j,i)=dft_ratio26{n}(5);%j=1,i=2,..8 order  j,i =5*8 for contourf
                DFT_vel{n}(j,i)=dft_ratio26_vel{n}(5);
                DFT_acc{n}(j,i)=dft_ratio26_acc{n}(5);
                p_mfr{n}(j,i)=dft_p26{n}(5);
                p_vel{n}(j,i)=dft_p26_vel{n}(5);
                p_acc{n}(j,i)=dft_p26_acc{n}(5);
%                 Hil{i,j}=Hilbert{5};
            elseif (j==5)
                DFT_mfr{n}(j,i)=dft_ratio26{n}(1);%j=5,i=2,..8
                DFT_vel{n}(j,i)=dft_ratio26_vel{n}(1);
                DFT_acc{n}(j,i)=dft_ratio26_acc{n}(1);
                p_mfr{n}(j,i)=dft_p26{n}(1);
                p_vel{n}(j,i)=dft_p26_vel{n}(1);
                p_acc{n}(j,i)=dft_p26_acc{n}(1);
%                 Hil{i,j}=Hilbert{1};
            else
                DFT_mfr{n}(j,i)=dft_ratio26{n}(intersect(find(plot_col==i),find(plot_row==j))); %j=2,3,4,i=2..8
                DFT_vel{n}(j,i)=dft_ratio26_vel{n}(intersect(find(plot_col==i),find(plot_row==j)));
                DFT_acc{n}(j,i)=dft_ratio26_acc{n}(intersect(find(plot_col==i),find(plot_row==j)));
                p_mfr{n}(j,i)=dft_p26{n}(intersect(find(plot_col==i),find(plot_row==j)));
                p_vel{n}(j,i)=dft_p26_vel{n}(intersect(find(plot_col==i),find(plot_row==j)));
                p_acc{n}(j,i)=dft_p26_acc{n}(intersect(find(plot_col==i),find(plot_row==j)));
%                 Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
            end
        else
            DFT_mfr{n}(j,i)=dft_ratio26{n}(intersect(find(plot_col==i),find(plot_row==j))); %i=1, j=1 2 3 4 5
            DFT_vel{n}(j,i)=dft_ratio26_vel{n}(intersect(find(plot_col==i),find(plot_row==j)));
            DFT_acc{n}(j,i)=dft_ratio26_acc{n}(intersect(find(plot_col==i),find(plot_row==j)));
            p_mfr{n}(j,i)=dft_p26{n}(intersect(find(plot_col==i),find(plot_row==j)));
            p_vel{n}(j,i)=dft_p26_vel{n}(intersect(find(plot_col==i),find(plot_row==j)));
            p_acc{n}(j,i)=dft_p26_acc{n}(intersect(find(plot_col==i),find(plot_row==j)));
%             Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:cellno
for i = 1:8
    for j= 1:5
        if p_mfr{n}(j,i) < 0.05
            p_mfr_p{n}(j,i)=1;
        else
            p_mfr_p{n}(j,i)=0;
        end
        if p_vel{n}(j,i) < 0.05
            p_vel_p{n}(j,i)=1;
        else
            p_vel_p{n}(j,i)=0;
        end
        if p_acc{n}(j,i) < 0.05
            p_acc_p{n}(j,i)=1;
        else
            p_acc_p{n}(j,i)=0;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change order 270-225-180-135-90-45-0-315
    for n=1:cellno
     for i=1:9;%length(unique_azimuth)+1                 %
        for j=1:5;%length(unique_elevation)
            if (i < 8 )                                 
               DFT_mfr270{n}(j,i)=DFT_mfr{n}(j,8-i); 
               DFT_vel270{n}(j,i)=DFT_vel{n}(j,8-i);
               DFT_acc270{n}(j,i)=DFT_acc{n}(j,8-i);
               p_mfr_p270{n}(j,i)=p_mfr_p{n}(j,8-i);
               p_vel_p270{n}(j,i)=p_vel_p{n}(j,8-i);
               p_acc_p270{n}(j,i)=p_acc_p{n}(j,8-i);
               
            elseif(i==8)
               DFT_mfr270{n}(j,i)=DFT_mfr{n}(j,i); 
               DFT_vel270{n}(j,i)=DFT_vel{n}(j,i);
               DFT_acc270{n}(j,i)=DFT_acc{n}(j,i);
               p_mfr_p270{n}(j,i)=p_mfr_p{n}(j,i);
               p_vel_p270{n}(j,i)=p_vel_p{n}(j,i);
               p_acc_p270{n}(j,i)=p_acc_p{n}(j,i);
%                 new_d_p_acc(i, j)=d_p_acc(i, j);
            else
               DFT_mfr270{n}(j,i)=DFT_mfr{n}(j,7); 
               DFT_vel270{n}(j,i)=DFT_vel{n}(j,7);
               DFT_acc270{n}(j,i)=DFT_acc{n}(j,7);
               p_mfr_p270{n}(j,i)=p_mfr_p{n}(j,7);
               p_vel_p270{n}(j,i)=p_vel_p{n}(j,7);
               p_acc_p270{n}(j,i)=p_acc_p{n}(j,7);
            end
        end
     end
    end
%%%%%   plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f=1:ceil(cellno./4);% ceil (46 cells/4) is 11.5 turned to 12 figures
    figure(f+5);orient landscape;
    for p=1:4;
        
            subplot(6,4,(p));contourf( DFT_mfr270{4*(f-1)+p}(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                 set(gca, 'xtick', [] );title( name2(4*(f-1)+p) );
                 set(gca, 'ytick', [] ); 
            subplot(6,4,(p+4));contourf( p_mfr_p270{4*(f-1)+p}(:,:) );colorbar('YTickLabel',{'NotSig',' ','Signf'});set(gca, 'ydir' , 'reverse');
                set(gca, 'xtick', [] );t=['pMFR=',num2str( P_anova(4*(f-1)+p) ),' / #Sig=',num2str( count_sig_dft(4*(f-1)+p) )];title(t);
                set(gca, 'ytick', [] ); 
            subplot(6,4,(p+8));contourf( DFT_vel270{4*(f-1)+p}(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                set(gca, 'xtick', [] );
                set(gca, 'ytick', [] ); 
            subplot(6,4,(p+12));contourf( p_vel_p270{4*(f-1)+p}(:,:) );colorbar('YTickLabel',{'NotSig',' ','Signf'});set(gca, 'ydir' , 'reverse');
                set(gca, 'xtick', [] );t=['pVel=',num2str( P_anova_vel(4*(f-1)+p) ),' / #Sig=',num2str( count_sig_dft_vel(4*(f-1)+p) )];title(t);
                set(gca, 'ytick', [] ); 
            subplot(6,4,(p+16));contourf( DFT_acc270{4*(f-1)+p}(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                set(gca, 'xtick', [] );
                set(gca, 'ytick', [] ); 
            subplot(6,4,(p+20));contourf( p_acc_p270{4*(f-1)+p}(:,:) );colorbar('YTickLabel',{'NotSig',' ','Signf'});set(gca, 'ydir' , 'reverse');
                set(gca, 'xtick', [] );t=['pAcc=',num2str( P_anova_acc(4*(f-1)+p) ),' / #Sig=',num2str( count_sig_dft_acc(4*(f-1)+p) )];title(t);
                set(gca, 'ytick', [] ); 
        
    end
end