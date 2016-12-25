%----------------------------------------------------------
%find variance/mean fit for d' stat
%for direction tuning - Z. Briggs
%---------------------------------

function var_mean_fit(batchfiledir, filename)
Curvefit_defines;
ploton=1;
printon=0;
pause on;

backdoor_fit_file= 'C:\MATLAB6p5\work\tempo_backdoor\results.txt';
backdoor_load_dir = 'C:\MATLAB6p5\work\tempo_backdoor';
data_output_file = 'C:\MATLAB6p5\work\tempo_backdoor\var_mean_mat.mat';
backdoor_load_ext = '_load.mat';
% 
% backdoor_fit_file= 'Z:\Users\Zack\tempo_backdoor\results.txt';
% backdoor_load_dir = 'Z:\Users\Zack\tempo_backdoor';
% data_output_file = 'Z:\Users\Zack\tempo_backdoor\var_mean_mat.mat';
% backdoor_load_ext = '_load.mat';


%Resolution of fit
% azstep=2*pi/63;
% elstep=pi/31;
azstep=2*pi/127;
elstep=pi/63;

el=-pi/2:elstep:pi/2;
az=-pi/2:azstep:3*pi/2;

if printon==1
    ploton=1;
end
%mapping for half rectified model
curvefit_parameter_mapping(CURVEFIT_PARAM_AZIMUTH) = 1;
curvefit_parameter_mapping(CURVEFIT_PARAM_ELEVATION) = 2;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN) = 3;
curvefit_parameter_mapping(CURVEFIT_PARAM_AMPLITUDE) = 4;
curvefit_parameter_mapping(CURVEFIT_PARAM_DC_OFFSET) = 5;
curvefit_parameter_mapping(CURVEFIT_PARAM_WEIGHT180) = 6;
curvefit_parameter_mapping(CURVEFIT_PARAM_NONLIN2) = 7;

errcnt=0;
cnt = 1;  %keeps track of cell number (not including those that are not skipped)
totalcnt=1; %keeps track of all cells including those that are skipped (for use with 
err=' ';


Dpopx=0;
Dpopy=0;
Tpopx=0;

filename = fullfile(batchfiledir, filename);
fid = fopen(filename);
% load('C:\matlab6p5\work\tempo_backdoor\r2chi2.mat');
line = fgetl(fid);

x_min = -90; x_max = 270; y_min = -90; y_max = 90;
x_tick = [-90:45:270]; y_tick = [-90:45:90];
%get bestfit params from results.txt

[file_names fitting_model_str fitting_bound_str tuning_model_str ...
        file_stims file_num_params residual x0str all_residual_str plot_r2_x plot_chi2_x] = textread( backdoor_fit_file, '%s %s %s %s %d %d %f %[^|]| %[^|]| %f %f' );
    
    num_files = length(file_names);
    
    while (line ~= -1)
        %pause
        % format for batch files
        % PATH  FILE    
        
        if (line(1) ~= '%')
            err=' ';
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
            backdoor_load_file = fullfile(backdoor_load_dir, [FILE backdoor_load_ext]);
            load(backdoor_load_file,'-mat');
            % go through the backdoor
            %             if cnt==1
            % convert the azimuth and elevation into a format amenable to plotting
            unique_azimuth = sort(unique_azimuth)' - 90;  % is the sort redundant?
            unique_elevation = sort(unique_elevation)';
            unique_azimuth = [unique_azimuth unique_azimuth(1) + 360];  % wrap around sphere
            [plot_azimuth plot_elevation] = meshgrid(unique_azimuth, unique_elevation);
            %             end
            %check for at least 3 repetitions
            skip=0;                

            %set a variable to grab the right gaze angle if there are more
            %than one.  AdDed by GCD, 1/21/05
            if (length(resp_mat(2,:,1,1)) > 1)  % a variable gaze run
                gazeindex = 2;  %use middle gaze angle
            else
                gazeindex = 1;
            end
            
            for m=1:26 
                if length(resp_mat(2,gazeindex,m,:))<4   %less then four because first value is the mean
                    skip=1;
                end
            end
            
            
            if skip==0
                
                for m=1:26
                    m_stat(cnt,m)=resp_mat(2,gazeindex,m,1);
                    L=length(resp_mat(2,gazeindex,m,:));
                    var_stat(cnt,m)=var(resp_mat(2,gazeindex,m,2:L));
                    if var_stat(cnt,m)==0
                        var_stat(cnt,m)=.1;
                        sprintf('%s %s', 'Replaced var=0 with var=.1 - Cell - ', FILE)
                        err=' Replaced var=0 with var=.1'
                        
                    end
                    
                    if m_stat(cnt,m)==0
                        m_stat(cnt,m)=.1;
                        sprintf('%s %s', 'Replaced mean=0 with mean=.1 - Cell -', FILE)
                        err= ' Replaced mean=0 with mean=.1'
                    end
                    lm(cnt,m)=log(m_stat(cnt,m));
                    lv(cnt,m)=log(var_stat(cnt,m));
                end
                
                %generate fit
                if ~(any(isinf(lm(cnt,:))) | any(isinf(lv(cnt,:))))
                    
                    [b(:,cnt),stat(:,cnt)]=robustfit(lm(cnt,:),lv(cnt,:));
                    %                                 stat(:,cnt)
                                        
                    %Calc diff
                    select = ( strcmp(FILE, file_names));
                    n = find(select == 1);
                    n = n(length(n));
                    num_params = file_num_params(n);
                    [best_x(1:num_params,cnt), num] = sscanf( x0str{n}, '%f', inf );
                    r2(cnt) = plot_r2_x(n);
                    chi2(cnt) = plot_chi2_x(n);               
                    if num < num_params
                        ['Error reading x0 parameters "' x0str{n} '"']
                        continue;
                    end
                    
                    for i=1:length(el)
                        for j=1:length(az)
                            F(i,j,cnt)= Curvefit_cos_tuning_7p_halfrectmodel_for_dprime(best_x(1:num_params,cnt),az(j),el(i));
                            
                        end
                    end
                    
                    for i=1:length(unique_elevation)
                        for j=1:length(unique_azimuth)
                            Fit_for_res(i,j,cnt) = Curvefit_cos_tuning_7p_halfrectmodel_for_dprime(best_x(1:num_params,cnt),(unique_azimuth(j)/180)*pi,(unique_elevation(i)/180)*pi);
                        end
                    end
                    %taken from tuning function
                    if ~all(F(:,:,cnt) == 0)
                        tmpmax=max(max(abs(F(:,:,cnt))));
                        F(:,:,cnt) = F(:,:,cnt) ./ tmpmax;
                        Fit_for_res(:,:,cnt) = Fit_for_res(:,:,cnt) ./ max(max(abs(Fit_for_res(:,:,cnt))));
                    end
                    F(:,:,cnt) = best_x(4,cnt).*F(:,:,cnt) + best_x(5,cnt);    %bestx 4 = amplitude
                    Fit_for_res(:,:,cnt) = best_x(4,cnt).*Fit_for_res(:,:,cnt) + best_x(5,cnt);    %bestx 4 = amplitude
                    
                    for i=1:length(F(:,1,cnt))
                        for j=1:length(F(1,:,cnt))
                            if(F(i,j,cnt)<0)
                                F(i,j,cnt)=0;
                            end
                        end
                    end
                    
                    for i=1:length(Fit_for_res(:,1,cnt))
                        for j=1:length(Fit_for_res(1,:,cnt))
                            if(Fit_for_res(i,j,cnt)<0)
                                Fit_for_res(i,j,cnt)=0;
                            end
                        end
                    end
                    
                    %log0 = inf .. problems?
                    var_fit(:,:,cnt) = exp(log(F(:,:,cnt)).*b(2,cnt) + b(1,cnt));
                    %                                                              for i=1:length(var_fit(:,1,cnt))
                    %                                                                  for j=1:length(var_fit(1,:,cnt))
                    %                                                                      if var_fit(i,j,cnt)==0
                    %                                                                          var_fit(i,j,cnt)=NaN;
                    %                                                                          sprintf('%s %s','Replaced var_fit=0 with var_fit=NaN cell-', FILE)
                    %                                                                     end
                    %                                                                 end
                    %                                                             end
                    %calc gradient

                    
                    
                    [FX(:,:,cnt),FY(:,:,cnt)]=gradient(F(:,:,cnt),azstep,elstep); 
                    
                    
                    %Matlab estimates first and last column/row
                    %derivatives, correction assuming cosine wraparound
                    
                    FX(:,1,cnt)=(F(:,2,cnt)-F(:,length(az)-1,cnt))./(2*azstep);
                    FX(:,length(az),cnt)=FX(:,1,cnt);
                    %what about Y?
                    
                    %calc dprime
                    std_fit(:,:,cnt)=sqrt(abs(var_fit(:,:,cnt)));
                    L=length(el);
                    %use 2:L-1 to eliminate poles from computation
                    
                    
                    if (~any(any(std_fit(2:L-1,:,cnt)<.0001)))
                        
                        Dx(:,:,cnt)=abs((FX(2:(L-1),:,cnt)./std_fit(2:L-1,:,cnt)));  %eliminates data at -90 and 90
                        Tx(:,:,cnt)=abs((FX(:,:,cnt)./std_fit(:,:,cnt)));
                        Dy(:,:,cnt)=abs((FY(:,:,cnt)./std_fit(:,:,cnt)));
                        
                        %For minimum threshold in x direction vs. azimuth
                        %plot
                        warning off MATLAB:divideByZero;
                        xThresh(:,:,cnt)=1./Tx(:,1:(length(az)-1),cnt);
                        warning on MATLAB:divideByZero;
                        minxThresh(cnt)=min(min(xThresh(:,:,cnt)));
                        %                         select=(xThresh(:,:,cnt)==minxThresh(cnt));
                        [i,j]=find((xThresh(:,:,cnt)==minxThresh(cnt)));
                        %                         if (sum(sum(select))==0 | sum(sum(select))>1)
                        if ~(length(i)==1) | ~(length(j)==1)
                            sprintf('%s %s %f %f', 'Error in cell, minxThresh count ',FILE, length(i),length(j))
                            errcnt=errcnt+1;
                            err=strcat(err,' Error in cell, minxThresh count');
                            minxThreshAZ(cnt)=NaN;
                            minxThreshEL(cnt)=NaN;
                            do=0;
                        else
                            minxThreshAZ(cnt)=az(j);
                            minxThreshEL(cnt)=el(i);
                            do=1;
                        end
                        
                        if do==1;
                            if ~(j<3 | j>(length(az)-3))
                                
                                tmp=min(min(xThresh(1:length(el),1:j-2,cnt)));
                                
                                tmp2=min(min(xThresh(1:length(el),j+2:length(az)-1,cnt)));
                                
                                
                                
                                tmp=min([tmp tmp2]);
                                [tmpi,tmpj]=find(tmp==xThresh(:,:,cnt));
                                if length(tmpi)>1
                                    sprintf('%s','error with new dbl x thresh')
                                    minxThreshAZ2(cnt)=NaN;
                                    minxThreshEL2(cnt)=NaN;
                                else
                                    minxThreshAZ2(cnt)=az(tmpj);
                                    minxThreshEL2(cnt)=el(tmpi);
                                    minxThresh2(cnt)=tmp;
                                end
                                
                            else
                                if j<3
                                    tmp=min(min(xThresh(1:length(el),4:length(az)-1,cnt)));
                                else
                                    tmp=min(min(xThresh(1:length(el),1:length(az)-4,cnt)));
                                end
                                minxThresh2(cnt)=tmp;
                                
                                [tmpi,tmpj]=find(tmp==xThresh(:,:,cnt));
                                minxThreshAZ2(cnt)=az(tmpj);
                                minxThreshEL2(cnt)=el(tmpi);
                            end
                        end
                        
                        
                        warning off MATLAB:divideByZero;
                        yThresh(:,:,cnt)=1./Dy(:,:,cnt);
                        warning on MATLAB:divideByZero;
                        minyThresh(cnt)=min(min(yThresh(:,:,cnt)));
                        %                         select=(xThresh(:,:,cnt)==minxThresh(cnt));
                        [i,j]=find((yThresh(:,:,cnt)==minyThresh(cnt)));
                        %                         if (sum(sum(select))==0 | sum(sum(select))>1)
                        if ~(length(i)==1) | ~(length(j)==1)
                            if((i(1)==i(2)) & (j(1)==1 | j(1)==length(az)) & (j(2)==1|j(2)==length(az))) %if it wraps around
                                minyThreshAZ(cnt)=az(j(1));
                                minyThreshEL(cnt)=el(i(1));
                            else
                                sprintf('%s %s %f %f', 'Error in cell, minyThresh count ',FILE, length(i),length(j))
                                errcnt=errcnt+1;
                                err=strcat(err,' Error in cell, minyThresh count');
                                do=0;
                            end
                        else
                            minyThreshAZ(cnt)=az(j);
                            minyThreshEL(cnt)=el(i);
                            do=1;
                        end
                        if do==1
                            if ~(j<3 | j>(length(az)-3))
                                if(i>2)
                                    tmp=min(yThresh(1:i-2,1:length(az),cnt));
                                    tmp2=min(yThresh(i-2:i,1:j-2,cnt));
                                    tmp3=min(yThresh(i-2:i,j+2:length(az),cnt));
                                else
                                    tmp=[];
                                    tmp2=[];
                                    tmp3=[];
                                end
                                if i<(length(el)-2)
                                    tmp4=min(yThresh(i+2:length(el),1:length(az),cnt));
                                    tmp5=min(yThresh(i:i+2,1:j-2,cnt));
                                    tmp6=min(yThresh(i:i+2,j+2:length(az),cnt));
                                else
                                    tmp4=[];
                                    tmp5=[];
                                    tmp6=[];
                                end
                                tmp=min([tmp tmp2 tmp3 tmp4 tmp5 tmp6]);
                                [tmpi,tmpj]=find(tmp==yThresh(:,:,cnt));
                                if length(tmpi)>1
                                    sprintf('%s','error with new dbl y thresh')
                                    minyThreshAZ2(cnt)=NaN;
                                    minyThreshEL2(cnt)=NaN;
                                else
                                    minyThreshAZ2(cnt)=az(tmpj);
                                    minyThreshEL2(cnt)=el(tmpi);
                                end
                                
                            else
                                if (j<3)
                                    range=4:length(az);
                                else
                                    range=1:length(az)-4;
                                end
                                
                                tmp=min(yThresh(1:i-2,1:length(az),cnt));
                                tmp2=min(yThresh(i+2:length(el),1:length(az),cnt));
                                tmp3=min(yThresh(i-2:i+2,range,cnt));
                                tmp=min([tmp tmp2 tmp3]);
                                [tmpi,tmpj]=find(tmp==yThresh(:,:,cnt));
                                if length(tmpi)>1
                                    sprintf('%s', 'More then one secondary min Y thresh found')
                                    minxThreshAZ2(cnt)=NaN;
                                    minxThreshEL2(cnt)=NaN;
                                else
                                    minyThreshAZ2(cnt)=az(tmpj);
                                    minyThreshEL2(cnt)=el(tmpi);
                                end
                            end
                        end
                        
                        
                        
                        
                        
                        %Width of peak
                        maxX(cnt)=max(max(F(:,:,cnt)));
                        [i,j]=find(maxX(cnt)==F(:,:,cnt));
                        if (length(i)==0 | length(i)>1 | length(j)>1)
                            sprintf('%s %s %f', 'Error in cell, maxX count ', FILE, sum(select))
                            errcnt=errcnt+1;
                            err=strcat(err, ' Error in cell, maxX count');
                            widthAZ(cnt)=NaN;
                            widthEL(cnt)=NaN;
                            
                        else
                            maxEL(cnt)=el(i);
                            maxAZ(cnt)=az(j);
                            
                            %                         tmp=F(i,:,cnt)-min(F(i,:,cnt));
                            tmp=F(i,:,cnt)-best_x(5,cnt);
                            fifty=tmp(j)/2;
                            fifty=fifty+best_x(5,cnt);
                            %                         fifty=fifty+min(F(i,:,cnt));
                            %                         tmp=find(abs(fifty-F(i,:,cnt))<.09);
                            if az(j)>pi/2
                                step=-1;
                            else
                                step=1;
                            end
                            
                            %                             loc=j;
                            %                             go=1;
                            %                             dblcheck=1;
                            %                             while go==1 & loc>1 & loc<length(az) & dblcheck==1
                            %                                 loc=loc+step;
                            %                                 go=abs(F(i,loc,cnt)-fifty)>abs(F(i,loc+step,cnt)-fifty);
                            %                                 dblcheck=F(i,loc,cnt)>F(i,loc+step,cnt);    
                            %                             end
                            
                            loc=10*(j+2*step);
                            
                            go=1;
                            dblcheck=1;
                            tmp=interp(F(i,:,cnt),10);
                            while go==1 & loc>1 & loc<length(tmp) & dblcheck==1
                                loc=loc+step;
                                go=abs(tmp(loc)-fifty)>abs(tmp(loc+step)-fifty);
                                dblcheck=tmp(loc)>tmp(loc+step);    
                            end
                            
                            
                            if dblcheck==0
                                sprintf('%s', 'At local minimum x')
                                err=strcat(err, ' At local minimum x');
                                widthAZ(cnt)=NaN;
                            else
                                if go==1
                                    sprintf('%s %f', 'Error finding az width, loc = ', loc)
                                    widthAZ(cnt)=NaN;
                                else                                                     
                                    %                                     loc=loc/10;
                                    %                                     loc=round(loc);
                                    loc=(-pi/2-azstep)+(loc/10)*azstep;
                                    widthAZ(cnt)=2*abs(az(j)-loc);
                                end
                            end
                            
                            tmp=F(:,j,cnt)-best_x(5,cnt);
                            fifty=tmp(i)/2;
                            fifty=fifty+best_x(5,cnt);
                            %                         fifty=fifty+min(F(i,:,cnt));
                            %                         tmp=find(abs(fifty-F(i,:,cnt))<.09);
                            if el(i)>0
                                step=-1;
                            else
                                step=1;
                            end
                            
                            %                             loc=j;
                            %                             go=1;
                            %                             dblcheck=1;
                            %                             while go==1 & loc>1 & loc<length(az) & dblcheck==1
                            %                                 loc=loc+step;
                            %                                 go=abs(F(i,loc,cnt)-fifty)>abs(F(i,loc+step,cnt)-fifty);
                            %                                 dblcheck=F(i,loc,cnt)>F(i,loc+step,cnt);    
                            %                             end
                            
                            loc=10*(i+2*step);
                            
                            go=1;
                            dblcheck=1;
                            tmp=interp(F(:,j,cnt),10);
                            while go==1 & loc>1 & loc<length(tmp) & dblcheck==1
                                loc=loc+step;
                                go=abs(tmp(loc)-fifty)>abs(tmp(loc+step)-fifty);
                                dblcheck=tmp(loc)>tmp(loc+step);    
                            end
                            
                            
                            if dblcheck==0
                                sprintf('%s', 'At local minimum y')
                                err=strcat(err, ' At local minimum y');
                                widthEL(cnt)=NaN;
                            else
                                if go==1
                                    sprintf('%s %f', 'Error finding el width, loc = ', loc)
                                    widthEL(cnt)=NaN;
                                else                                                     
                                    %                                   loc=loc/10;
                                    %                                   loc=round(loc);
                                    loc=(-pi/2-elstep)+(loc/10)*elstep;
                                    widthEL(cnt)=2*abs(el(i)-loc);
                                end
                            end
                        end
                        
                        Dpopx=Dpopx+Dx(:,:,cnt).^2;
                        Tpopx=Tpopx+Tx(:,:,cnt).^2;
                        
                        %                         Dx(Dx<.0001)=NaN;
                        %                         dSx(:,:,cnt)=1./Dx(:,:,cnt);
                        
                        Dpopy=Dpopy+Dy(:,:,cnt).^2;
                        %                         Dy(Dy<.0001)=NaN;
                        %                         dSy(:,:,cnt)= 1./Dy;
                        plotdx=1;
                    else
                        plotdx=0;
                    end
                    
                    
                    name{cnt}=sprintf(FILE);            
                    
                    if ploton==1
                        %Do plots
                        
                        g=figure;
                        set(g, 'position', [10,10,1000,500])
                        subplot(2,3,1);
                        scatter(lm(cnt,:),lv(cnt,:));
                        hold on;
                        var_line(cnt,:) = lm(cnt,:)*b(2,cnt) + b(1,cnt);
                        plot(lm(cnt,:),var_line(cnt,:),'-');
                        xlabel('log(mean)');
                        ylabel('log(var)');
                        hold off;
                        
                        
                        %plot quiver/contour
                        subplot(2,3,3);
                        %                        contourf(F(:,:,cnt));
                        contourf(Fit_for_res(:,:,cnt));
                        
                        %Auto scale axis based on mean data while
                        %keeping the min value at zero
                        AX=caxis;
                        if AX(1)>0
                            AX(1)=0;
                        end
                        caxis(AX);
                        
                        set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                        xlabel('azimuth'); ylabel('elevation');
                        title('Fit Data');
                        colorbar;
                        %                     hold on;
                        %                     quiver(FX(:,:,cnt),FY(:,:,cnt));
                        %                     hold off;
                        %close(g);
                        for m=1:5
                            for n=1:9
                                % i can not think of a cleaner way to do this stupid
                                % all azimuths at the pole thing.
                                e = 1e-10;  % floating point tolerance
                                select = ((abs(unique_point_azimuth - mod(unique_azimuth(n),360)) < e | ...
                                    unique_elevation(m) == 90 | unique_elevation(m) == -90) & ...
                                    abs(unique_point_elevation - unique_elevation(m)) < e);
                                if sum(select) == 0;
                                    eat_shit_and_die = 1;
                                end
                                
                                plot_mean_data(m,n) = resp_mat(2,gazeindex,select,1);
                            end
                        end
                        Res(:,:,cnt)=Fit_for_res(:,:,cnt)-plot_mean_data(:,:);
                        
                        
                        % plot the mean data
                        
                        subplot(2,3,2);
                        %subplot('Position',[0 0 0.5 0.5]);
                        [C h] = contourf(plot_azimuth,plot_elevation,plot_mean_data);
                        caxis(AX);
                        colorbar;
                        set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                        axis([x_min,x_max,y_min,y_max]); 
                        xlabel('azimuth'); ylabel('elevation');
                        title('Mean Data');
                        
                        %                         title('Mean Data'); set( h,'Interpreter', 'none'); set(h,'FontName','FixedWidth'); .
                        if plotdx==1
                            subplot(2,3,4);
                            %quiver(Tx(:,:,cnt),Dy(:,:,cnt));
                            contourf(Dx(:,:,cnt));
                            set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                            xlabel('azimuth'); ylabel('elevation');
                            title('Dx');
                            colorbar;
                        end
                        subplot(2,3,5);
                        axis off;
                        text(0,1,strcat(FILE,err));
                        text(0,1-.1,sprintf('%s %f','Prefered AZimuth = ',(best_x(1,cnt))*180/pi));
                        text(0,1-.2,sprintf('%s %f','Prefered ELevation = ',(best_x(2,cnt)*180/pi)));
                        text(0,1-.3,sprintf('%s %f','Amplitude = ',best_x(4,cnt)));
                        text(0,1-.4,sprintf('%s %f','DC Offset = ',best_x(5,cnt)));
                        text(0,1-.5,sprintf('%s %f','Weight 180 = ',best_x(6,cnt)));
                        text(0,1-.6,sprintf('%s %f','Nonlin1 = ',best_x(3,cnt)));
                        text(0,1-.7,sprintf('%s %f','Nonlin2 = ',best_x(7,cnt)));
                        text(0,1-.8,sprintf('%s %f','R^2 value = ',r2(cnt)));
                        text(0,1-.9,sprintf('%s %f','Chi^2 value = ',chi2(cnt)));
                        
                        subplot(2,3,6);
                        %subplot('Position',[0 0 0.5 0.5]);
                        [C h] = contourf(plot_azimuth,plot_elevation,Res(:,:,cnt));
                        caxis(AX);
                        colorbar;
                        set(gca,'ydir','reverse','xdir','reverse','xtick',x_tick,'ytick',y_tick);
                        axis([x_min,x_max,y_min,y_max]); 
                        xlabel('azimuth'); ylabel('elevation');
                        title('Residual')
                        
                        
                        
                        pause;
                        if printon == 1
                            set(gcf,'PaperOrientation','Landscape', 'PaperPositionMode', 'Auto');
                            print -dwinc;
                        end
                        close(g);
                        
                    end
                    
                    cnt = cnt + 1; 
                    
                else
                    sprintf('%s %s %s', 'Skipped cell' , FILE, ' log of zero')
                    err=strcat(err, ' skipped, log of zero');
                end
                
            else
                sprintf('%s %s', 'Skipped cell, less then 3 reps - ', FILE)   
                err=strcat(err, ' skipped,less then 3 reps');
            end 
            
            errindex{totalcnt,1}=FILE;
            errindex{totalcnt,2}=err;
            totalcnt=totalcnt+1;
            
        end  
        
        line = fgetl(fid);
    end
    Dpopx=sqrt(Dpopx);
    Dpopy=sqrt(Dpopy);
       
    save(data_output_file, 'm_stat', 'var_stat','b','stat','F', 'best_x', 'FX', 'FY', 'var_fit', 'name', 'std_fit', 'Dpopx', 'Dpopy', 'Dx', 'Dy', 'lm', 'lv','Tpopx','Tx', 'minxThresh','minxThreshAZ', 'errindex', 'minyThresh','minyThreshAZ', 'minyThreshEL', 'minxThreshEL','widthAZ', 'widthEL','r2','chi2', 'minxThreshAZ2','minxThreshEL2', 'minxThresh2', 'minyThreshAZ2','minyThreshEL2','maxEL', 'maxAZ','az','el');
    
    
