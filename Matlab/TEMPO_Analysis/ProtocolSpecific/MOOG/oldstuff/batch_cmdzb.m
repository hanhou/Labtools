function batch_cmd(batchfiledir, filename)
pause off;
if nargin < 1 | strcmp(batchfiledir, '');
    batchfiledir = 'Z:\data\Tempo\Batch Files';
end
if nargin < 2 | strcmp(filename, '');
    filename = 'MOOG_vary_fix.m';
end

% do different global plotting stuff.
Curvefit_defines;

filename = fullfile(batchfiledir, filename);
fid = fopen(filename);

% copied from BATCH_GUI_Tempo_Analysis
line = fgetl(fid);
cnt = 0;
while (line ~= -1)
    
    %pause
    % format for batch files
    % PATH  FILE    
    
    if (line(1) ~= '%')
        
        % first remove any comment text at the end of a line (following a %), GCD, added 9/25/01
        comment_start = find(line == '%');
        if ~isempty(comment_start)
            line = line(1:(comment_start(1)-1));
        end
        
        spaces = isspace(line);
        space_index = find(spaces);
        
        %get path / file
        PATH = line(1:space_index(1) - 1);
        FILE = line(space_index(1) + 1:space_index(2) - 1)
        l = length(FILE);
        if (FILE(l-3:l) == '.htb')	% .htb extension already there
            filename = [PATH FILE];   %the HTB data file
            logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
        else	%no extension in FILE, add extensions
            filename = [PATH FILE '.htb'];   %the HTB data file
            logfile = [PATH FILE '.log'];   %the TEMPO log file
        end

        % go through the backdoor
        cnt = cnt + 1;
        Heading_CurveFit([],[],[],[],[],[],[],[],[],[],PATH,FILE);

%         % close the plots
%         %closehandle = figure;
%         for closeindex = [2 4 6]
%             close(closeindex);
%         end

        pause;

        select = (curvefit_fitting_stims_export == 2);
        if any(select)
                all_files{cnt} = FILE;
                plot_r2_x(cnt) = curvefit_R2_export(select,1);
                plot_chi2_x(cnt) = curvefit_Chi2_export(select,1);
%               plot_r2_y(cnt) = curvefit_R2_export(select,2);
% %             plot_r2_x2(cnt) = curvefit_R2_export(select,3);
%                 plot_r2_y2(cnt) = curvefit_R2_export(select,4);

%                  fic=fopen('C:\matlab6p5\work\r2.txt', 'a');   % open output file
%                  fprintf( fic, '\n%s %d %d', FILE, plot_r2_x(cnt), plot_chi2_x(cnt) );
%                  fprintf( fic, '\n' );
%                  fclose(fic);

            end
%             plot_r2_y3(cnt) = curvefit_R2_export(select,5);
%             plot_r2_x3(cnt) = curvefit_R2_export(select,6);

%             plot_az_y_horz(cnt,:) = [0 0 0];
%             plot_az_y_vert(cnt,:) = [0 0 0];
%             plot_el_y_horz(cnt,:) = [0 0 0];
%             plot_el_y_vert(cnt,:) = [0 0 0];
% 
%             gas=[-20 0 20];
%             for n=1:length(gas);
%                 selectg = (curvefit_gaze_angle_export == gas(n));
%                 m = find(selectg == 1);
%                 if any(selectg)
%                     if curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_X
%                         plot_az_y_horz(cnt,n) = curvefit_fitting_param_export(select,3,7*(m-1)+1);
%                         plot_el_y_horz(cnt,n) = curvefit_fitting_param_export(select,3,7*(m-1)+2);
%                     elseif curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_Y
%                         plot_az_y_vert(cnt,n) = curvefit_fitting_param_export(select,3,7*(m-1)+1);
%                         plot_el_y_vert(cnt,n) = curvefit_fitting_param_export(select,3,7*(m-1)+2);
%                     end
%                 end
%             end
        end
        
    
    line = fgetl(fid);
end
% save('C:\matlab6p5\work\tempo_backdoor\r2chi2.mat','plot_r2_x', 'plot_chi2_x');
% % figure(776);
% scatter(plot_r2_x, plot_r2_y);
% x=[0:0.1:1]; y=[0:0.1:1];
% hold on; plot(x,y,'r-', 'LineWidth', 2); hold off; 
% ylabel('R^{2} eye centered common pos gain');
% xlabel('R^{2} head centered common pos gain');
% 
% % figure(555);
% % scatter(plot_r2_x2, plot_r2_y2);
% % x=[0:0.1:1]; y=[0:0.1:1];
% % hold on; plot(x,y,'r-', 'LineWidth', 2); hold off; 
% % ylabel('R^{2} eye centered 5p common gain');
% % xlabel('R^{2} head centered 5p common gain');
% 
% diffp = plot_r2_x - plot_r2_y;
% %select = (diffp <= 0);
% %diffp(select) = 0;
% diffp
% plot_r2_x
% all_files

% figure(556);
% scatter(plot_r2_x3, plot_r2_y3);
% x=[0:0.1:1]; y=[0:0.1:1];
% hold on; plot(x,y,'r-', 'LineWidth', 2); hold off; 
% ylabel('R^{2} eye centered');
% xlabel('R^{2} head centered');
% 
% 
% % move az rotation angles to range [-pi/2 3*pi/2)
% select = plot_az_y_horz < -pi/2;
% plot_az_y_horz(select) = plot_az_y_horz(select) + 2*pi;
% select = plot_az_y_horz >= 3*pi/2;
% plot_az_y_horz(select) = plot_az_y_horz(select) - 2*pi;
% 
% select = plot_az_y_vert < -pi/2;
% plot_az_y_vert(select) = plot_az_y_vert(select) + 2*pi;
% select = plot_az_y_vert >= 3*pi/2;
% plot_az_y_vert(select) = plot_az_y_vert(select) - 2*pi;
% 
% % move el rotation angles to range [-pi/2 pi/2]
% select = plot_el_y_horz < -pi/2;
% plot_el_y_horz(select) = plot_el_y_horz(select) + pi;
% select = plot_el_y_horz < -pi/2;
% plot_el_y_horz(select) = plot_el_y_horz(select) + pi;
% select = plot_el_y_horz > pi/2;
% plot_el_y_horz(select) = plot_el_y_horz(select) - pi;
% select = plot_el_y_horz > pi/2;
% plot_el_y_horz(select) = plot_el_y_horz(select) - pi;
% 
% select = plot_el_y_vert < -pi/2;
% plot_el_y_vert(select) = plot_el_y_vert(select) + pi;
% select = plot_el_y_vert < -pi/2;
% plot_el_y_vert(select) = plot_el_y_vert(select) + pi;
% select = plot_el_y_vert > pi/2;
% plot_el_y_vert(select) = plot_el_y_vert(select) - pi;
% select = plot_el_y_vert > pi/2;
% plot_el_y_vert(select) = plot_el_y_vert(select) - pi;
% 
% % remove incomplete data sets
% a = [0 0 0 0];
% for n=1:cnt
%     if all(plot_az_y_horz(n,:) ~= 0)
%         a(1) = a(1) + 1;
%         plot_az_y_horz2(n,:) = plot_az_y_horz(n,:);
%     end
%     if all(plot_az_y_vert(n,:) ~= 0)
%         a(2) = a(2) + 1;
%         plot_az_y_vert2(n,:) = plot_az_y_vert(n,:);
%     end
%     if all(plot_el_y_horz(n,:) ~= 0)
%         a(3) = a(3) + 1;
%         plot_el_y_horz2(n,:) = plot_el_y_horz(n,:);
%     end
%     if all(plot_el_y_vert(n,:) ~= 0)
%         a(4) = a(4) + 1;
%         plot_el_y_vert2(n,:) = plot_el_y_vert(n,:);
%     end
% end
% 
% plot_az_y_horz2 = plot_az_y_horz2/pi*180;
% plot_az_y_vert2 = plot_az_y_vert2/pi*180;
% plot_el_y_horz2 = plot_el_y_horz2/pi*180;
% plot_el_y_vert2 = plot_el_y_vert2/pi*180;
% 
% diff(plot_az_y_horz2,2);
% 
% figure(777);
% subplot(2,2,1);
% plot(gas,plot_az_y_horz2);
% title('az rot versus gaze for horz fix');
% subplot(2,2,2);
% plot(gas,plot_az_y_vert2);
% title('az rot versus gaze for vert fix');
% subplot(2,2,3);
% plot(gas,plot_el_y_horz2);
% title('el rot versus gaze for horz fix');
% subplot(2,2,4);
% plot(gas,plot_el_y_vert2);
% title('el rot versus gaze for vert fix');
