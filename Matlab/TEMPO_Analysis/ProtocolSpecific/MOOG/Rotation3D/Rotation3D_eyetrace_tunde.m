function Rotation3D_eyetrace_tunde(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for figure writing , if running batch, comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % From here RVOR/pursuit
% if size(data.eye_data,1)>6
%     LEFT_EYE_1_2=9;
%     RIGHT_EYE_3_4=10;
% else
%     LEFT_EYE_1_2=7;
%     RIGHT_EYE_3_4=8;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %************************************************************************%
% %plot Vertical vs. Horizontal
% switch (data.eye_flag)
%     case (LEFT_EYE_1_2)
%         Eye_Select='Left Eye'
%         Hor=1;        Ver=2;
%     case(RIGHT_EYE_3_4)
%         Eye_Select='Right Eye'
%         Hor=3;        Ver=4;
% end
% %

% % set process type
% process_type = 'tempo_gui';
process_type = 'local_mat_files';
% % % % % % % % % % % % % % % % % % 

if strcmp(process_type, 'local_mat_files' ) == 1
    FILE_mat = FILE;
    ext_position = strfind(FILE, '.');
    FILE_mat( (ext_position + 1) : end ) = 'mat'; %change extention to 'mat'    
    save(FILE_mat)
    disp('matlab file saved.')
else %were running through tempo_gui ...


    % Pursuit_eachCell_eyetrace



    temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
    temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
    temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
    temp_spike_rates = data.spike_rates(SpikeChan, :);
    temp_total_trials = data.misc_params(OUTCOME, :);
    temp_spike_data = data.spike_data(1,:);   % spike rasters
    % temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);

    % null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
    null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );
    %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
    % trials = 1:length(temp_azimuth);
    trials = 1:length(temp_elevation);
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    azimuth = temp_azimuth(~null_trials & select_trials);
    elevation = temp_elevation(~null_trials & select_trials);
    stim_type = temp_stim_type(~null_trials & select_trials);
    spike_rates = temp_spike_rates(~null_trials & select_trials);
    % fp_rotate = temp_fp_rotate(~null_trials & select_trials);

    unique_azimuth = munique(azimuth');
    unique_elevation = munique(elevation');
    unique_stim_type = munique(stim_type');
    % unique_fp_rotate = munique(fp_rotate');

    % repeat = floor( length(temp_spike_rates) / (length(unique_stim_type)*(length(unique_fp_rotate)*(length(unique_elevation)-2)+2)+1) );
    %Yong's bad Aihua's good for calculate repeat

    trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_stim_type) + 1;
    repetition = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

    %%%%%%%%%%%%%%  clear   %%%%%%%%%%%%%%%%%%%%%%%%%
    clear offset_x offset_y resp_x resp_y resp_x_up resp_y_up resp_x_down resp_y_down resp_x_left resp_y_left resp_x_right resp_y_right
    clear eye_x_up eye_y_up eye_x_down eye_y_down eye_x_left eye_y_left eye_x_right eye_y_right mean_eye_x_up mean_eye_y_up mean_eye_x_down mean_eye_y_down mean_eye_x_left mean_eye_y_left mean_eye_x_right mean_eye_y_right
    clear vis_eye_x_up vis_eye_y_up vis_eye_x_down vis_eye_y_down vis_eye_x_left vis_eye_y_left vis_eye_x_right vis_eye_y_right vis_mean_eye_x_up vis_mean_eye_y_up vis_mean_eye_x_down vis_mean_eye_y_down vis_mean_eye_x_left vis_mean_eye_y_left vis_mean_eye_x_right vis_mean_eye_y_right

    for k = 1 : length(unique_stim_type)
        {'k =', k, 'unique_stim_type = ' , unique_stim_type(k)}
        for i = 1 : length(unique_azimuth)
            {'i =', i, 'unique_azi_type =', unique_azimuth(i)}
            for j = 1 : length(unique_elevation)
                {'j =', j, 'unique_ele', unique_elevation(j)}
                select = find( temp_azimuth==unique_azimuth(i) & temp_elevation==unique_elevation(j) & temp_stim_type==unique_stim_type(k) );
                %           the above line "select = find ..." basically finds the indices/index where all the conditions are met
                %           i.e. the index/indices where temp_azimuth is equal to unique_azimuth(i), and so on... ---- Tunde
                if sum(select)>0
                    for jj = 1 : repetition % for convenience with sacrefice of some of the trials
                        {'jj =', jj, 'repetition'}
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%        select eye Right ch.3, 4 or Left ch.1, 2
                        %%%%
                        % % %                Right Eye Old_Que, Zeblon, Lothar, Que_Laby, Zebu_Laby
                        %                     offset_x = mean( data.eye_data(3,201:300,select(jj)) ); % horizontal
                        %                     offset_y = mean( data.eye_data(4,201:300,select(jj)) ); % vertical
                        %                     resp_x{k,i,j}(jj,:) = data.eye_data(3,201:600,select(jj)) - offset_x;  % horizontal
                        %                     resp_y{k,i,j}(jj,:) = data.eye_data(4,201:600,select(jj)) - offset_y;  %

                        %                Left Eye Azrael, Lothar(after..), New_Que, but some eye data are destroyed !! use old cell for Que (see Translation)
                        %% and Lother left Eye
                        offset_x = mean( data.eye_data(1,201:300,select(jj)) ); % horizontal
                        offset_y = mean( data.eye_data(2,201:300,select(jj)) ); % vertical
                        resp_x{k,i,j}(jj,:) = data.eye_data(1,201:600,select(jj)) - offset_x;  % horizontal; 400pts == 2 secs
                        resp_y{k,i,j}(jj,:) = data.eye_data(2,201:600,select(jj)) - offset_y;  % vertical; 400pts == 2 secs


                        %                Eye select (for figure wrinting)
                        %                     offset_x = mean( data.eye_data(Hor,201:300,select(jj)) ); % horizontal
                        %                     offset_y = mean( data.eye_data(Ver,201:300,select(jj)) ); % vertical
                        %                     resp_x{k,i,j}(jj,:) = data.eye_data(Hor,201:600,select(jj)) - offset_x;  % horizontal
                        %                     resp_y{k,i,j}(jj,:) =
                        %                     data.eye_data(Ver,201:600,select(jj)) - offset_y;  %
                        %                     vertical
                    end
                else
                    resp_x{k,i,j}(:,:) = resp_x{k,1,j}(:,:);%??? +-90 elevation
                    resp_y{k,i,j}(:,:) = resp_y{k,1,j}(:,:);
                    %                 this results in the matrices for  pages 1 and/or 5 for
                    %                 the 3D cell array, having the same values. See scanned
                    %                 notes...Tunde
                end
            end
        end
        resp_x_up{k}(:,:) = resp_x{k,1,3}(:,:);     resp_y_up{k}(:,:) = resp_y{k,1,3}(:,:);
        resp_x_down{k}(:,:) = resp_x{k,5,3}(:,:);   resp_y_down{k}(:,:) = resp_y{k,5,3}(:,:);
        resp_x_left{k}(:,:) = resp_x{k,1,1}(:,:);   resp_y_left{k}(:,:) = resp_y{k,1,1}(:,:);
        resp_x_right{k}(:,:) = resp_x{k,1,5}(:,:);  resp_y_right{k}(:,:) = resp_y{k,1,5}(:,:);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Difference between mean for repeatition and 201-600 should be
    %%%%%%%% attention!!
    % paper=0;

    for k = 1 : length(unique_stim_type)

        eye_x_up(:,:)= resp_x_up{k}(:,101:300);% vestibular (k=1) only
        eye_y_up(:,:) = resp_y_up{k}(:,101:300);%(jj, 201-600) jj=repeatition, 400 points(1point=50ms)
        eye_x_down(:,:)= resp_x_down{k}(:,101:300);%middle 1 sec = 101:300
        eye_y_down(:,:) = resp_y_down{k}(:,101:300);
        eye_x_left(:,:)= resp_x_left{k}(:,101:300);
        eye_y_left(:,:)= resp_y_left{k}(:,101:300);
        eye_x_right(:,:)= resp_x_right{k}(:,101:300);
        eye_y_right(:,:)= resp_y_right{k}(:,101:300);

        mean_eye_x_up{k} = mean(eye_x_up(:,:)');
        mean_eye_y_up{k} = mean(eye_y_up(:,:)');
        mean_eye_x_down{k} = mean(eye_x_down(:,:)');
        mean_eye_y_down{k} = mean(eye_y_down(:,:)');
        mean_eye_x_left{k} = mean(eye_x_left(:,:)');
        mean_eye_y_left{k} = mean(eye_y_left(:,:)');
        mean_eye_x_right{k} = mean(eye_x_right(:,:)');
        mean_eye_y_right{k} = mean(eye_y_right(:,:)');

        %%%%% here plot figures %%%% commented
        % figure(k+10);
        % subplot(2,2,1);
        % plot(mean_eye_x_up{k} ,'rx');hold on;plot(mean_eye_x_down{k} ,'ro');hold off;title('Hor / pitch Up (x) and Down(o)');ylim([-10 10]);
        % subplot(2,2,2);
        % plot(mean_eye_x_left{k} ,'rx');hold on;plot(mean_eye_x_right{k} ,'ro');hold off;title('Hor / yaw Left (x) and Right(o)');ylim([-10 10]);
        % subplot(2,2,3);
        % plot(mean_eye_y_up{k} ,'bx');hold on;plot(mean_eye_y_down{k} ,'bo');hold off;title('Ver / pitch Up (x) and Down(o)');ylim([-10 10]);
        % subplot(2,2,4);
        % plot(mean_eye_y_left{k} ,'bx');hold on;plot(mean_eye_y_right{k} ,'bo');hold off;title('Ver / yaw Left (x) and Right(o)');ylim([-10 10]);

    end


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%   take mean to plotsupplemental figure, use folder 'analysis' mat
    % %%%%%%   code to calculate
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear re_x_up re_y_up re_x_down re_y_down re_x_left re_y_left re_x_right re_y_right

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     mean for 5 repetition for figures

    for k = 1 : length(unique_stim_type)
        re_x_up{k} = mean(resp_x_up{k}(:,:));
        re_y_up{k} = mean(resp_y_up{k}(:,:));
        re_x_down{k} = mean(resp_x_down{k}(:,:));
        re_y_down{k} = mean(resp_y_down{k}(:,:));
        re_x_left{k} = mean(resp_x_left{k}(:,:));
        re_y_left{k} = mean(resp_y_left{k}(:,:));
        re_x_right{k} = mean(resp_x_right{k}(:,:));
        re_y_right{k} = mean(resp_y_right{k}(:,:));

        %     the followig code is done to scale the figures
        eyescale{k}=[re_x_up{k} re_y_up{k} re_x_down{k} re_y_down{k} re_x_left{k} re_y_left{k} re_x_right{k} re_y_right{k}];
        eyemax{k}=max(eyescale{k})
        eyemin{k}=min(eyescale{k})
        if eyemax{k}<0.5
            eyemax{k}=0.5;
        end
        if eyemin{k}>-0.5
            eyemin{k}=-0.5;
        end
    end
    %


    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %For figure writing, in order to batch comment!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


    title1 = 'up';
    title2 = 'down';
    title3 = 'left';
    title4 = 'right';

    subtitle{1}='Vestibular';subtitle{2}='Visual';

    lengthstim=length(unique_stim_type);
    if lengthstim > 2; % combined No!
        lengthstim=2;
    end

    %%%%%%%%%%%%%% MEAN Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tunde1 = figure(3);  orient landscape;
    axes('position',[0 0 1 1]);
    xlim([1,10]);
    ylim([1,10]);
    axis off
    %     for k = 1 : length(unique_stim_type)
    for k = 1 : lengthstim;% combined No!
        text (3+5*(k-1), 9, subtitle(k));
    end
    text (5, 9, FILE);

    %     for k = 1 : length(unique_stim_type)
    for k = 1 : lengthstim
        for i = 1:4

            axes('position',[0.1+0.5*(k-1) 0.65-0.2*(i-1) 0.35 0.15]);

            if i==1
                plot(re_x_up{k},'r.');
                hold on;
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');

                ylabel('(deg)');
                title(['Eye Position /  ',title1]);

                plot(re_y_up{k},'b.');
                hold off;

            elseif i==2
                plot(re_x_down{k},'r.');
                hold on;
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title2]);

                plot(re_y_down{k},'b.');
                hold off;

            elseif i==3
                plot(re_x_left{k},'r.');
                hold on;
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title3]);

                plot(re_y_left{k},'b.');
                hold off;

            elseif i==4
                plot(re_x_right{k},'r.');
                hold on;
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title4]);

                plot(re_y_right{k},'b.');
                hold off;


            end
        end
    end

    %%%%%%%%%%%%%% Each repetition Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tunde2 = figure(4);  orient landscape;
    axes('position',[0 0 1 1]);
    xlim([1,10]);
    ylim([1,10]);
    axis off
    %     for k = 1 : length(unique_stim_type)
    for k = 1 : lengthstim;% combined No!
        text (3+5*(k-1), 9, subtitle(k));
    end
    text (5, 9, FILE);
    %  %%%%%%%%%%%%%%%%%%%%% max_ and min_ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %
    %  still working on it
    %
    %
    %     for k = 1 : length(unique_stim_type)
    %
    %     max_re_x_up{k} = max(max(resp_x_up{k}(:,:)));
    %     max_re_y_up{k} = max(max(resp_y_up{k}(:,:)));
    %     max_re_x_down{k} = max(max(resp_x_down{k}(:,:)));
    %     max_re_y_down{k} = max(max(resp_y_down{k}(:,:)));
    %     max_re_x_left{k} = max(max(resp_x_left{k}(:,:)));
    %     max_re_y_left{k} = max(max(resp_y_left{k}(:,:)));
    %     max_re_x_right{k} = max(max(resp_x_right{k}(:,:)));
    %     max_re_y_right{k} = max(max(resp_y_right{k}(:,:)));
    %
    %     max_eyescale{k}=[max_re_x_up{k} max_re_y_up{k} max_re_x_down{k} max_re_y_down{k} max_re_x_left{k} max_re_y_left{k} max_re_x_right{k} max_re_y_right{k}];
    %     max_eyemax{k}=max(max_eyescale{k})
    %
    %     min_eyescale{k}=[min_e_x_up{k} min_re_y_up{k} min_re_x_down{k} min_re_y_down{k} min_re_x_left{k} min_re_y_left{k} min_re_x_right{k} min_re_y_right{k}];
    %     min_eyemin{k}=min(min_eyescale{k})
    %
    %     if eyemax{k}<0.5
    %         eyemax{k}=0.5;
    %     end
    %     if eyemin{k}>-0.5
    %         eyemin{k}=-0.5;
    %     end
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    howmany=size(resp_x_up{1})

    %     for k = 1 : length(unique_stim_type)
    for k = 1 : lengthstim; % permit only 2 conditions, because of lack of space
        for i = 1:4

            axes('position',[0.1+0.5*(k-1) 0.65-0.2*(i-1) 0.35 0.15]);

            if i==1
                for jj=1:howmany(1)
                    plot(resp_x_up{k}(jj,:),'r.');
                    hold on;
                    plot(resp_y_up{k}(jj,:),'b.');
                    hold on;
                end
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');

                ylabel('(deg)');
                title(['Eye Position /  ',title1]);
                %                  for jj=1:howmany(1)
                %                  plot(resp_y_up{k}(jj,:),'b.');
                %                     hold on;
                %                 end

            elseif i==2
                for jj=1:howmany(1)
                    plot(resp_x_down{k}(jj,:),'r.');
                    hold on;
                    plot(resp_y_down{k}(jj,:),'b.');
                    hold on;
                end
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title2]);

                %                 plot(resp_y_down{k},'b.');
                %                     hold off;

            elseif i==3
                for jj=1:howmany(1)
                    plot(resp_x_left{k}(jj,:),'r.');
                    hold on;
                    plot(resp_y_left{k}(jj,:),'b.');
                    hold on;
                end
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title3]);

                %                     plot(resp_y_left{k},'b.');
                %                          hold off;

            elseif i==4
                for jj=1:howmany(1)
                    plot(resp_x_right{k}(jj,:),'r.');
                    hold on;
                    plot(resp_y_right{k}(jj,:),'b.');
                    hold on;
                end
                xlim( [1, 400] );
                ylim( [eyemin{k}, eyemax{k}] );
                set(gca, 'XTickMode','manual');
                set(gca, 'xtick',[1,100,200,300,400]);
                set(gca, 'xticklabel','0|0.5|1|1.5|2');
                ylabel('(deg)');
                title(['Eye Position /  ',title4]);

                %                     plot(resp_y_right{k},'b.');
                %                         hold off;


            end
        end
    end



    saveas(tunde1, [FILE, '1.fig'], 'fig')



end

% for k = 1 : length(unique_stim_type)
%
%     figure(k+20)
%
% subplot(4,1,1)
% plot(re_x_up{k},'r.');
%     hold on;
%     xlim( [1, 400] );
%     ylim( [eyemin{k}, eyemax{k}] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2');
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title1]);
%
%     plot(re_y_up{k},'b.');
%     hold off;
%
%
% subplot(4,1,2)
% plot(re_x_down{k},'r.');
%     hold on;
%     xlim( [1, 400] );
%     ylim( [eyemin{k}, eyemax{k}] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2');
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title2]);
%
%     plot(re_y_down{k},'b.');
%     hold off;
%
% subplot(4,1,3)
% plot(re_x_left{k},'r.');
%     hold on;
%     xlim( [1, 400] );
%     ylim( [eyemin{k}, eyemax{k}] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2');
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title3]);
%
%     plot(re_y_left{k},'b.');
%     hold off;
%
% subplot(4,1,4)
% plot(re_x_right{k},'r.');
%     hold on;
%     xlim( [1, 400] );
%     ylim( [eyemin{k}, eyemax{k}] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2');
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title4]);
%
%     plot(re_y_right{k},'b.');
%     hold off;
%
% end
%
%
%
% %     %% commented 04/04/07
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%                         output to text file
% % sprint_txt = ['%s'];
% % for i = 1 : 400*repeat*10
% %     sprint_txt = [sprint_txt, ' %1.2f'];
% % end
% %
% %
% % % if you want to save 'world', select 2, 'pursuit', select 3
% %
% %
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
% %                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% % % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
% % %                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% % %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% % %
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\QueEye_Rot_Vet.dat'];
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\AzraelEye_Rot_Vet.dat'];
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\AzraelEye_Rot_Vis.dat'];
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\ZebulonEye_Rot_Vet.dat'];
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\QueLaby_Rot_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\ZebulonLaby_Rot_Vet.dat'];
% %
% % printflag = 0;
% % if (exist(outfile, 'file') == 0)    %file does not yet exist
% %     printflag = 1;
% % end
% % fid = fopen(outfile, 'a');
% % if (printflag)
% %     fprintf(fid, 'FILE\t');
% %     fprintf(fid, '\r\n');
% % end
% % fprintf(fid, '%s', buff);
% % fprintf(fid, '\r\n');
% % fclose(fid);
% %
% %
% %
% % %%%%%                         output to text file
% % sprint_txt2 = ['%s'];
% % for i = 1 : 400*repeat*10
% %     sprint_txt2 = [sprint_txt2, ' %1.2f'];
% % end
% %
% %
% % % if you want to save 'world', select 2, 'pursuit', select 3
% %
% %
% % % buff= sprintf(sprint_txt2, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
% % %                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
% %                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% % %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% % %
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\AzraelEye_Rot_Vet.dat'];
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\AzraelEye_Rot_Vis.dat'];
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\QueEye_Rot_Vis.dat'];
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\ZebulonEye_Rot_Vis.dat'];
% %
% % % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\QueLaby_Rot_Vis.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\ZebulonLaby_Rot_Vis.dat'];
% %
% % printflag = 0;
% % if (exist(outfile, 'file') == 0)    %file does not yet exist
% %     printflag = 1;
% % end
% % fid = fopen(outfile, 'a');
% % if (printflag)
% %     fprintf(fid, 'FILE\t');
% %     fprintf(fid, '\r\n');
% % end
% % fprintf(fid, '%s', buff);
% % fprintf(fid, '\r\n');
% % fclose(fid);



return;

