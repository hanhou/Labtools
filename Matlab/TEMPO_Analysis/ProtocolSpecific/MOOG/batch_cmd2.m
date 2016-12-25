function batch_cmd(batchfiledir, filename)

if nargin < 1 | strcmp(batchfiledir, '');
    batchfiledir = 'Z:\data\Tempo\Batch Files';
end
if nargin < 2 | strcmp(filename, '');
    filename = 'MOOG_vary_fix.m';
end

% do different global plotting stuff.
Curvefit_defines;

% some inits
plot_eq_x = [0];
plot_eq_y = [0];
plot_horz_sighead_x = [0];
plot_horz_sighead_y = [0];
plot_horz_sigeye_x = [0];
plot_horz_sigeye_y = [0];
plot_horz_nsighead_x = [0];
plot_horz_nsighead_y = [0];
plot_horz_nsigeye_x = [0];
plot_horz_nsigeye_y = [0];
plot_vert_sighead_x = [0];
plot_vert_sighead_y = [0];
plot_vert_sigeye_x = [0];
plot_vert_sigeye_y = [0];
plot_vert_nsighead_x = [0];
plot_vert_nsighead_y = [0];
plot_vert_nsigeye_x = [0];
plot_vert_nsigeye_y = [0];

filename = fullfile(batchfiledir, filename);
fid = fopen(filename);

% copied from BATCH_GUI_Tempo_Analysis
line = fgetl(fid);
cnt = 0;
good_cnt = 0;
good_cnt_horz = 0;
good_cnt_vert = 0;
gas_vert = zeros(1,3);
gas_horz = zeros(1,3);
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
        OPT_FILE = fullfile('C:\MATLAB6p5\work\tempo_backdoor', 'Curvefit_default_options.m');

        % go through the backdoor
        cnt = cnt + 1;
        DirectionTuningPlot_Curvefit([],[],[],[],[],[],[],[],[],[],PATH,FILE,OPT_FILE);

%         % close the plots
%         %closehandle = figure;
%         for closeindex = [2 4 6]
%             close(closeindex);
%         end

        %pause;

        select = (curvefit_fitting_stims_export == 2);
        if any(select)
            good_cnt = good_cnt + 1;
            all_files_x{cnt} = FILE;
            
            if curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_X
                good_cnt_horz = good_cnt_horz+1;
                plot_horz_x(cnt) = curvefit_R2_export(select,1);
                plot_horz_y(cnt) = curvefit_R2_export(select,2);
%                 a_horz(good_cnt_horz,1:3) = [curvefit_fitting_param_export(select,1,5) ...
%                         curvefit_fitting_param_export(select,1,8) curvefit_fitting_param_export(select,1,11)];
%                 d_horz(good_cnt_horz,1:3) = [curvefit_fitting_param_export(select,1,6) ...
%                         curvefit_fitting_param_export(select,1,9) curvefit_fitting_param_export(select,1,12)];
%                 b_horz(good_cnt_horz,1:3) = [curvefit_fitting_param_export(select,1,7) ...
%                         curvefit_fitting_param_export(select,1,10) curvefit_fitting_param_export(select,1,13)];
%                 gas_horz(good_cnt_horz,1:3) = 10000*ones(1,3);
%                 gas_horz(good_cnt_horz,1:length(curvefit_gaze_angle_export)) = curvefit_gaze_angle_export';
%                 ratio_param_horz(good_cnt_horz) = curvefit_fitting_param_export(select,3,5);
                plot_horz_x_diff(cnt) = curvefit_R2_export(select,2) - curvefit_R2_export(select,1);
%                 plot_horz_k(good_cnt_horz) = curvefit_fitting_param_export(select,3,5);
            elseif curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_Y
                good_cnt_vert = good_cnt_vert+1;
                plot_vert_x(cnt) = curvefit_R2_export(select,1);
                plot_vert_y(cnt) = curvefit_R2_export(select,2);
%                 a_vert(good_cnt_vert,1:3) = [curvefit_fitting_param_export(select,1,5) ...
%                         curvefit_fitting_param_export(select,1,8) curvefit_fitting_param_export(select,1,11)];
%                 d_vert(good_cnt_vert,1:3) = [curvefit_fitting_param_export(select,1,6) ...
%                         curvefit_fitting_param_export(select,1,9) curvefit_fitting_param_export(select,1,12)];
%                 b_vert(good_cnt_vert,1:3) = [curvefit_fitting_param_export(select,1,7) ...
%                         curvefit_fitting_param_export(select,1,10) curvefit_fitting_param_export(select,1,13)];
%                 gas_vert(good_cnt_vert,1:3) = 10000*ones(1,3);
%                 gas_vert(good_cnt_vert,1:length(curvefit_gaze_angle_export)) = curvefit_gaze_angle_export';
%                 ratio_param_vert(good_cnt_vert) = curvefit_fitting_param_export(select,3,5);
                plot_vert_x_diff(cnt) = curvefit_R2_export(select,2) - curvefit_R2_export(select,1);
%                 plot_vert_k(good_cnt_vert) = curvefit_fitting_param_export(select,3,5);
            end
            
%             if curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_EQUAL_MEANS
%                 plot_eq_x(cnt+1) = curvefit_R2_export(select,1);
%                 plot_eq_y(cnt+1) = curvefit_R2_export(select,2);
%             elseif curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_X
%                 if curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_1_SIG
%                     plot_horz_sighead_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_horz_sighead_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_2_SIG
%                     plot_horz_sigeye_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_horz_sigeye_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_1_NSIG
%                     plot_horz_nsighead_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_horz_nsighead_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_2_NSIG
%                     plot_horz_nsigeye_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_horz_nsigeye_y(cnt+1) = curvefit_R2_export(select,2);
%                 end
%             elseif curvefit_fixation_type_export == CURVEFIT_VARY_FIXATION_Y
%                 if curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_1_SIG
%                     plot_vert_sighead_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_vert_sighead_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_2_SIG
%                     plot_vert_sigeye_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_vert_sigeye_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_1_NSIG
%                     plot_vert_nsighead_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_vert_nsighead_y(cnt+1) = curvefit_R2_export(select,2);
%                 elseif curvefit_R2_distrib_sig_export(select,1) == CURVEFIT_R2_DISTRIB_2_NSIG
%                     plot_vert_nsigeye_x(cnt+1) = curvefit_R2_export(select,1);
%                     plot_vert_nsigeye_y(cnt+1) = curvefit_R2_export(select,2);
%                 end
%             end
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

% figure(776);
% x=[0:0.1:1]; y=[0:0.1:1];
% plot(x,y,'r-', 'LineWidth', 2);
% 
% if length(plot_eq_x) > 0
%     hold on; h0 = scatter(plot_eq_x, plot_eq_x,'x','g'); hold off;
% end
% plot_eq_str = sprintf('%d eq means', sum(plot_eq_x > 0));
% 
% if length(plot_horz_sighead_x) > 0
%     hold on; h1 = scatter(plot_horz_sighead_x, plot_horz_sighead_y, '+', 'b'); hold off;
% end
% plot_horz_sighead_str = sprintf('%d horz head sig', sum(plot_horz_sighead_x > 0));
% 
% if length(plot_horz_sigeye_x) > 0
%     hold on; h2 = scatter(plot_horz_sigeye_x, plot_horz_sigeye_y, 'o', 'b'); hold off;
% end
% plot_horz_sigeye_str = sprintf('%d horz eye sig', sum(plot_horz_sigeye_x > 0));
% 
% if length(plot_horz_nsighead_x) > 0
%     hold on; h3 = scatter(plot_horz_nsighead_x, plot_horz_nsighead_y, '+', 'r'); hold off;
% end
% plot_horz_nsighead_str = sprintf('%d horz head nsig', sum(plot_horz_nsighead_x > 0));
% 
% if length(plot_horz_nsigeye_x) > 0
%     hold on; h4 = scatter(plot_horz_nsigeye_x, plot_horz_nsigeye_y, 'o', 'r'); hold off;
% end
% plot_horz_nsigeye_str = sprintf('%d horz eye nsig', sum(plot_horz_nsigeye_x > 0));
% 
% if length(plot_vert_sighead_x) > 0
%     hold on; h5 = scatter(plot_vert_sighead_x, plot_vert_sighead_y, 's', 'b'); hold off;
% end
% plot_vert_sighead_str = sprintf('%d vert head sig', sum(plot_vert_sighead_x > 0));
% 
% if length(plot_vert_sigeye_x) > 0
%     hold on; h6 = scatter(plot_vert_sigeye_x, plot_vert_sigeye_y, '*', 'b'); hold off;
% end
% plot_vert_sigeye_str = sprintf('%d vert eye sig', sum(plot_vert_sigeye_x > 0));
% 
% if length(plot_vert_nsighead_x) > 0
%     hold on; h7 = scatter(plot_vert_nsighead_x, plot_vert_nsighead_y, 's', 'r'); hold off;
% end
% plot_vert_nsighead_str = sprintf('%d vert head nsig', sum(plot_vert_nsighead_x > 0));
% 
% if length(plot_vert_nsigeye_x) > 0
%     hold on; h8 = scatter(plot_vert_nsigeye_x, plot_vert_nsigeye_y, '*', 'r'); hold off;
% end
% plot_vert_nsigeye_str = sprintf('%d vert eye sig', sum(plot_vert_nsigeye_x > 0));
% 
% ylabel('R^{2} eye centered');
% xlabel('R^{2} head centered');
% legend([h0(1);h1(1);h2(1);h3(1);h4(1);h5(1);h6(1);h7(1);h8(1)], ...
%     plot_eq_str,plot_horz_sighead_str,plot_horz_sigeye_str,plot_horz_nsighead_str, ...
%     plot_horz_nsigeye_str,plot_vert_sighead_str,plot_vert_sigeye_str,...
%     plot_vert_nsighead_str,plot_vert_nsigeye_str,2);
% title(sprintf('cell total = %d', good_cnt));

% figure(234);
% subplot(1,2,1);
% hist(ratio_param_horz,10);
% ylabel('cnt');
% xlabel('k gaze angle ratio');
% title('horz');
% subplot(1,2,2);
% hist(ratio_param_vert,10);
% ylabel('cnt');
% xlabel('k gaze angle ratio');
% title('vert');
% 
% figure(235);
% subplot(1,2,1);
% scatter(plot_horz_x_diff, plot_horz_k);
% xlabel('R^{2} eye - R{2} head');
% ylabel('k gaze angle ratio');
% title('horz');
% subplot(1,2,2);
% scatter(plot_vert_x_diff, plot_vert_k);
% xlabel('R^{2} eye - R{2} head');
% ylabel('k gaze angle ratio');
% title('vert');

figure(236);
x=[0:0.1:1]; y=[0:0.1:1];
plot(x,y,'r-', 'LineWidth', 2);
hold on; scatter(plot_horz_x, plot_horz_y); hold off;
hold on; scatter(plot_vert_x, plot_vert_y, 's', 'r'); hold off;
ylabel('R^{2} eye centered gf');
xlabel('R^{2} eye centered');

plot_vert_x_diff
plot_horz_x_diff
all_files_x
find(plot_horz_x_diff == min(plot_horz_x_diff))
find(plot_vert_x_diff == min(plot_vert_x_diff))

% figure(237); clf;
% subplot(1,3,1);
% plot(gas_horz',a_horz', 'bx-');
% hold on; plot(gas_vert',a_vert','rx-');
% set(gca, 'XLim', [-30 30]);
% subplot(1,3,2);
% plot(gas_horz',d_horz', 'bx-');
% hold on; plot(gas_vert',d_vert','rx-');
% set(gca, 'XLim', [-30 30]);
% subplot(1,3,3);
% plot(gas_horz',b_horz','bx-');
% hold on; plot(gas_vert',b_vert','rx-');
% set(gca, 'XLim', [-30 30]);

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
