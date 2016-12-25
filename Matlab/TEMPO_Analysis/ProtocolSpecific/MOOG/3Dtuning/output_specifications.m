% % %%%%%%%% Output_specifications file***************
% % Please use this to independently save your output files to different locations based the the person running the script in question.
% % If you have any questions, please see Tunde

% Basic example
% % if (getenv('computername'), 'COMPUTERNAME')
% %     outfile = ..........
% % end

% % % % % % % % % % % % % % % % % % % % % % % % % % %AYANNA'S SECTION % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if (strcmp(getenv('computername'), 'AYANNA'))
    
    if (strcmp(Analysis{1}, 'Plot Rotation Tuning 3D'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI, var_term);
        outfile = ['Z:\Users\Ayanna\tempo_output\3D_Rotation_Tuning_Ayanna.dat'];

        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t');
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    
    %%%%% next %%%%%
    elseif (strcmp(Analysis{1}, 'Plot Tuning Surface'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std,  F_val, P_anova, DDI, var_term );

        outfile = ['Z:\Users\Ayanna\tempo_output\3D_Direction_Tuning_Ayanna.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.

        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t Veb_DDI\t Vis_DDI\t Com_DDI\t Veb_var_term\t Vis_var_term\t Com_var_term\t');
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t'); % gain variable removed!!! ab 29/9/07
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
        
    %%%%%next%%%%%%%    
    elseif (strcmp(Analysis{1}, 'Output Rotation PSTH(Katsu)'))

        buff = sprintf(sprint_txt, FILE, pickupmin{:} );
        outfile = ['Z:\Users\Ayanna\tempo_output\fake_rotation_PSTH_output_Ayanna2.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.

        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);

    %%%%%next%%%%    
    elseif  (strcmp(Analysis{1}, 'Output Translation PSTH(Katsu)'))

        buff = sprintf(sprint_txt, FILE, pickupmin{:} );
        outfile = ['Z:\Users\Ayanna\tempo_output\translation_PSTH_output_Ayanna2.dat']; %saves data to my z:\users folder ab 21 Sept. 2007.

        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % TUNDE'S SECTION% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if (strcmp(getenv('computername'), 'DATAPRIME1'))
    if (strcmp(Analysis{1}, 'Output Rotation PSTH (Katsu)'))      
        buff = sprintf(sprint_txt, FILE, pickupmax{:} );
        outfile = ['Z:\Users\Tunde\tempo_output\fake_rotation_PSTH_output_Ayanna.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.
        
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % %SYED's SECTION % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if (strcmp(getenv('computername'), 'SYED'))

    if (strcmp(Analysis{1}, 'Plot Rotation Tuning 3D'))

        buff = sprintf('%s\t %4.2f\t  %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t   %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI);
        % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\rotation3d\syed_Rotation_Tuning_Sum.dat'];
%         outfile=['F:\Syed_Analysis\syed_Rotation_Tuning_Sum.dat']
        outfile = ['Z:\Users\Syed Chowdhury\syed2\tempo_output\3D_Direction_Tuning_Syed.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Veb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Veb_HTI\t Veb_HTIerr\t Veb_P\t Veb_std\t gain\t Veb_F_anova\t Veb_anova\t Veb_DDI');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);


    end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % %HUI's SECTION % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
if (strcmp(getenv('computername'), 'HUI'))

    if (strcmp(Analysis{1}, 'Plot Tuning Surface'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std,  F_val, P_anova, DDI, var_term );

        outfile = ['Z:\Users\HuiM\tempo_output\3D_Direction_Tuning.dat']; %saves data to my z:\users folder ab 20 Sept. 2007.

        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            %     fprintf(fid, 'FILE\t         SPon\t Veb_min\t Vis_min\t Comb_min\t Veb_max\t Vis_max\t Comb_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Comb_azi\t Comb_ele\t Comb_amp\t Veb_HTI\t Vis_HTI\t Comb_HTI\t Veb_HTIerr\t Vis_HTIerr\t Comb_HTIerr\t Veb_P\t Vis_P\t Comb_P\t Veb_std\t Vis_std\t Comb_std\t gain\t F_anova\t P_anova\t Veb_DDI\t Vis_DDI\t Com_DDI\t Veb_var_term\t Vis_var_term\t Com_var_term\t');
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t'); % gain variable removed!!! ab 29/9/07
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    elseif (strcmp(Analysis{1}, 'Plot Rotation Tuning 3D'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI, var_term);
        outfile = ['Z:\Users\HuiM\tempo_output\3D_Rotation_Tuning.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);        
    elseif (strcmp(Analysis{1}, 'Output Firing rate'))
        buff = sprintf(sprint_txt, FILE, spon_resp, ves, ves_std);        
        outfile = ['Z:\Users\HuiM\tempo_output\Translation_FiringRate.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)   % file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);                 
    elseif(strcmp(Analysis{1}, 'Output ROT_Firing rate'))
        buff = sprintf(sprint_txt, FILE, ves, ves_std );
%         buff = sprintf(sprint_txt, FILE, ves );
        outfile = ['Z:\Users\HuiM\tempo_output\VesOnly_Rotation_FiringRate.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)   % file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);        
    end    
    
end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %Aihua's SECTION % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % %
if (strcmp(getenv('computername'), 'AIHUA'))

    if (strcmp(Analysis{1}, 'Plot Tuning Surface'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std,  F_val, P_anova, DDI, var_term,P_anova_hor);
        outfile = ['C:\Aihua\z_TempOutputs\3D_Direction_Tuning_Aihua.dat']; 
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t Veb_p_hor\t Vis_p_hor\t'); % gain variable removed!!! ab 29/9/07
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
%         saveas(2,['C:\Aihua\z_TempOutputs\figures\' FILE(1:end-4)],'png');               
        %saveas(2,['C:\Aihua\z_TempOutputs\figures\' FILE(1:end-4) '_Tuning'],'png');               
    elseif (strcmp(Analysis{1}, 'Plot Rotation Tuning 3D'))

        buff = sprintf('%s\t %4.2f\t   %4.3f\t   %4.3f\t   %4.3f\t   %4.3f\t  %4.3f\t  %4.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %6.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %1.3f\t  %1.3f\t  %1.3f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t  %2.4f\t', ...
            FILE, spon_resp, Min_resp, Max_resp, Vec_sum{:}, r, HTI_boot, p{:} , resp_std, gain, F_val, P_anova, DDI, var_term);
        outfile = ['C:\Aihua\z_TempOutputs\3D_Rotation_Tuning_Aihua.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t SPon\t Veb_min\t Vis_min\t Veb_max\t Vis_max\t Veb_azi\t Veb_ele\t Veb_amp\t Vis_azi\t Vis_ele\t Vis_amp\t Veb_HTI\t Vis_HTI\t Veb_HTIerr\t Vis_HTIerr\t Veb_P\t Vis_P\t Veb_std\t Vis_std\t gain\t Veb_F_anova\t Vis_F_anova\t Veb_P_anova\t Vis_P_anova\t Veb_DDI\t Vis_DDI\t Veb_var_term\t Vis_var_term\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
        saveas(2,['C:\Aihua\z_TempOutputs\figures\' FILE(1:end-4)],'png'); 
    elseif (strcmp(Analysis{1}, 'Output Firing rate'))
         buff = sprintf(sprint_txt, FILE, spon_resp, ves, ves_std); 
        outfile = ['C:\Aihua\z_TempOutputs\VesOnly_Translation_FiringRate.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)   % file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);        
    elseif(strcmp(Analysis{1}, 'Output ROT_Firing rate'))
        buff = sprintf(sprint_txt, FILE, ves, ves_std );
%         buff = sprintf(sprint_txt, FILE, ves );
        outfile = ['C:\Aihua\z_TempOutputs\VesOnly_Rotation_FiringRate.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)   % file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    end
end
