function RFContourFigurePlot(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;
Path_Defs;

%*************************************************************************%
g_spikeTracker=zeros(30,30,30,360);
g_dirCounter=zeros(30,30,360);
g_spikeTracker_f=zeros(30,30,30,360);
g_dirCounter_f=zeros(30,30,360);

%Grab the number of width and height patches
numXpatches=data.revcorr_params(PATCH_DIMS,1,1);
numYpatches=data.revcorr_params(PATCH_DIMS,1,2);

%Grab some parameters of the stimulus
View_Dist=data.one_time_params(VIEW_DIST);%View Distance
Screen_Xsiz=data.one_time_params(SCREEN_XSIZ);%Screen size: x
Screen_Ysiz=data.one_time_params(SCREEN_YSIZ);%Screen size: y
FieldLoc_x=0;% FieldLoc_x=unique(data.moog_params(FIELDLOC_X,:,1));%field location: x
FieldLoc_y=0;% FieldLoc_y=unique(data.moog_params(FIELDLOC_Y,:,1));%field location: y
FieldSize_x=Screen_Xsiz;% FieldSize_x=unique(data.moog_params(FIELDSIZE_X,:,1));%field size: x
FieldSize_y=Screen_Ysiz;% FieldSize_y=unique(data.moog_params(FIELDSIZE_Y,:,1));%field size: y

%Determine total number of directions per patch
temp_directions = squeeze(data.revcorr_params(DOT_DIREC, :, :));
x=find(isnan(temp_directions(1,:))==0);
directions=temp_directions(:,x);
uniqueDirs=munique(directions(1,:)')';
numDirections=length(directions(1,:))/(numXpatches*numYpatches);

for (Xpatch=1:numXpatches) 
    for (Ypatch=1:numYpatches)
        clear temp;temp=((Xpatch-1)*numDirections*numYpatches+(Ypatch-1)*numDirections);
        directions_array{Xpatch,Ypatch}=directions(:,temp+1:temp+numDirections);        
    end
end

%Set the correlation delay parameters
corrDelayLow=0;
corrDelayInc=1;
corrDelayHigh=200;

%Store all the different spike histograms we collect in this data structure
for (trial=1:size(data.revcorr_params,2))
    %Create the correlation delay matrix 
    corrDelay = [corrDelayLow : corrDelayInc : corrDelayHigh];
    
    %Pull out the original spikes
    origSpikes = squeeze(data.spike_data(SpikeChan,:,trial));  
    
    %figure out what the edges will be for the spike histogram
    histEdges = find(squeeze(data.spike_data(2,:,trial))~=0);
    meanEdge = round(mean(diff(histEdges)));
    histEdges = [histEdges,histEdges(length(histEdges)) + meanEdge];
    
    spikeHistIndex=0;
    %Do a spike histogram for all the correlation delays
    for (delay = corrDelay)
        spikeHistIndex = spikeHistIndex+1;
        %Shift the spikes over to compensate for neuron lag. Pad the end
        %with zeros because shifting will truncate the vector
        shiftedSpikes = origSpikes(delay+1:length(origSpikes));
        shiftedSpikes = [shiftedSpikes, zeros(1,5000-length(shiftedSpikes))];
        
        %Find all shiftedSpikes elements that aren't  zero
        spikeTimes=find(shiftedSpikes);
        
        %Do a histogram to bin out spike counts for each directions
        %presented in the trial
        tmpHist=histc(spikeTimes,histEdges);  tmpHist=tmpHist./[diff(histEdges) meanEdge]*1000;
        spikeHist{spikeHistIndex}(trial,:)=tmpHist(1:length(tmpHist)-1);
        
        %forward to get the mean noise level
        shiftedSpikes_f=zeros(1,delay+1);
        shiftedSpikes_f = [shiftedSpikes_f,origSpikes(1:5000-length(shiftedSpikes_f))];       
        spikeTimes_f= find(shiftedSpikes_f);
        tmpHist_f = histc(spikeTimes_f, histEdges);tmpHist_f=tmpHist_f./[diff(histEdges) meanEdge]*1000;
        spikeHist_f{spikeHistIndex}(trial,:) = tmpHist_f(1:length(tmpHist_f)-1);           
    end    
end

%OK. The basic strategy here is to iterate through every patch and use the
%histogram we did above to add on directions to make another histrogram,
%which will be our analysis.
bootstp_num=0;%bootstp_num=1000;%Do the bootstrap
[g_spikeTracker, maxDelay]=DirRFplot(spikeHist,directions_array,numXpatches,numYpatches,spikeHistIndex,uniqueDirs, corrDelay,PATH,FILE,bootstp_num,2);
[g_spikeTracker_f, maxDelay_f]=DirRFplot(spikeHist_f,directions_array,numXpatches,numYpatches,spikeHistIndex,uniqueDirs, corrDelay,PATH,FILE,bootstp_num,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%At max delay, fit direction vs response data for each loc
clear DelayIndex;DelayIndex=find(corrDelay==maxDelay);
%cplot_wrapped=zeros(numXpatches+2, numYpatches+2);
cplot_wrapped=zeros(numXpatches, numYpatches);
param_label = 'A          mu         sigma       DC';
MaxCount_Fitted=0;
figure(4);clf;orient landscape;figure(5);clf;orient landscape;
for Xpatch=1:numXpatches
    for Ypatch=1:numYpatches
        clear Response;Response=squeeze(g_spikeTracker(Xpatch,Ypatch,DelayIndex,uniqueDirs+1))';
        k=Xpatch + (Ypatch-1)*numXpatches;
        [Response_F, x, x0]=OrienFitting(Response);
        NewX{k}=x;
        NewX0{k}=x0;
        %xDirection=[0    45    90   135   180   225   270   315 ]*pi/180;
        xDirection=[0    45    90   135   180   225   270   315   360]*pi/180; Response=[Response Response(1)]; Response_F=[Response_F Response_F(1)];
        figure(4);subplot(numXpatches,numYpatches, Xpatch + (Ypatch-1)*numXpatches);plot(xDirection,Response,'bo');hold on; % ylim([0 maxCount])        
        temp_Ori1(k,:)=xDirection;
        temp_Ori2(k,:)=Response;        
        % finish plotting
        x_smooth = 0:0.01:2*pi;
        y_smooth = (DirCurvefit(x,x_smooth));
        xdata_tran = [0    45    90   135   180   225   270   315   360]*pi/180;%xdata_tran = [0    45    90   135   180   225   270   315]*pi/180;%xdata_tran = [-90 -45 0 45 90 135 180 225 270] * pi/180;
        x_smooth_tran = xdata_tran(1):0.01:xdata_tran(end);
        y_smooth_tran = (DirCurvefit(x,x_smooth_tran));
        plot(x_smooth_tran,y_smooth_tran,'r');        
        temp_Fit1(k,:)=x_smooth_tran;
        temp_Fit2(k,:)=y_smooth_tran;
        % take peak of fitted function as preferred heading 
        peak = x_smooth(find(y_smooth == max(y_smooth))) * 180/pi;
        cplot_wrapped(Ypatch,Xpatch)=max(y_smooth_tran)-min(y_smooth_tran);%cplot_wrapped(Xpatch,Ypatch)=max(y_smooth_tran)-min(y_smooth_tran);cplot_wrapped
        tempMaxCount_Fitted=max(y_smooth_tran);
        if (tempMaxCount_Fitted > MaxCount_Fitted)
            MaxCount_Fitted=tempMaxCount_Fitted
        end     
        xlim( [xdata_tran(1) xdata_tran(end)]);
        set(gca, 'xtick', xdata_tran);
        set(gca, 'xticklabel','0|45|90|135|180|225|270|315|360');  
        %set(gca, 'xticklabel','0|45|90|135|180|225|270|315');  
        Fit_Error(Xpatch,Ypatch)=100/length(xDirection)*sum(abs(Response_F-Response)/max(Response));        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get the circular variance at the maxDelay
        clear Rate Resp;
        Resp=sum(Response.*exp(2*i*xDirection))/sum(Response);
        CV_Max(Xpatch,Ypatch)=1-abs(Resp);   
        
        TempResp(Xpatch,Ypatch)=abs(Resp);
        [Azi_Max(Xpatch,Ypatch), Ele_Max(Xpatch,Ypatch), Amp_Max(Xpatch,Ypatch)] = vectorsum(Response/max(Response));
        r_Max(Xpatch,Ypatch) = 1-abs(Amp_Max(Xpatch,Ypatch)) / sum(Response/max(Response)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %check the goodness of fitting data
        clear FitData;FitData=DirCurvefit(x,xDirection);
        figure(5);subplot(numXpatches,numYpatches, Xpatch + (Ypatch-1)*numXpatches);plot(Response,FitData,'*');axis([0 200 0 200]);

        %Linear Regression
        [p1 mu]=polyfit(Response,FitData,1);
%        x2 = 0:.1*max(xDirection):max(xDirection);  
        x2=0:.1*max(XLim):max(XLim);
        y2=polyval(p1,x2); hold on;plot(x2,y2,'g'); 
        [R,p2] = corrcoef(Response,FitData);
        text(0.5*max(XLim),0.95*max(YLim),['Y= ' num2str(p1(1)) '*X + ' num2str(p1(2),'%4.3f')])
        text(0.5*max(XLim),0.85*max(YLim),['r= ' num2str(R(1,2),'%4.3f')])
        text(0.5*max(XLim),0.75*max(YLim),['p= ' num2str(p2(1,2),'%4.3f')])
        text(max(x2),max(y2),' \leftarrow Fitted Line','FontSize',8)
        R_CurveFit(Xpatch,Ypatch)=R(1,2);
      
        %Also, write out some summary data to a cumulative summary file
        sprint_txt=['%s\t %f\t   %f\t   %f\t   %f\t   %f\t   %f\t   %f\t   %f\t'];
        buff=sprintf(sprint_txt, FILE, x0(1), x0(2)*180/pi, x0(3)*180/pi, x0(4), x(1),x(2)*180/pi, x(3)*180/pi, x(4));
        outfile =['C:\Aihua\z_TempOutputs\CurveFit.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t         A0\t mu0\t sigma0\t DC0\t A\t mu\t sigma\t DC\t');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
        fclose(fid);
    end
end 

for (Xpatch = 1:(numXpatches*numYpatches))   
    figure(4);subplot(numXpatches, numYpatches, Xpatch);  ylim([0 MaxCount_Fitted+0.2]);
    % show params for each fit
    param_text = [num2str(NewX{Xpatch}(1),'%1.2f') '     ' num2str(NewX{Xpatch}(2)*180/pi,'%1.2f') '    ' num2str(NewX{Xpatch}(3)*180/pi,'%1.2f') '     ' num2str(NewX{Xpatch}(4),'%1.2f')];            
    x0_text = [num2str(NewX0{Xpatch}(1),'%1.2f') '     ' num2str(NewX0{Xpatch}(2)*180/pi,'%1.2f') '    ' num2str(NewX0{Xpatch}(3)*180/pi,'%1.2f') '     ' num2str(NewX0{Xpatch}(4),'%1.2f')];    
    y_lim = ylim;
    y_range = y_lim(2)-y_lim(1);
    text(0, y_lim(2)+.17*y_range, param_label,'FontSize', 8);
    text(0, y_lim(2)+.10*y_range, x0_text, 'FontSize', 8);
    text(0, y_lim(2)+.03*y_range, param_text, 'FontSize',8);     
    if (Xpatch == numXpatches*numYpatches)
        legend('Original Data', 'Fitted Curve',4);xlabel('Azimuth');ylabel('Normalized Response');       
    elseif (Xpatch == 1)
        % Puts a label that says the filename above all the subplots.
        text(-2, y_lim(2)+.35*y_range, [PATH, FILE]); text(750,4,['Time= -',num2str(maxDelay),'ms']);       
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check significance of Circular Variance and calculate p value, do bootstrap at the same
%time to test value varience
% bootstp_num=1000;
% [g_spikeTracker_f, maxDelay_f]=DirRFplot(spikeHist_f,directions_array,numXpatches,numYpatches,spikeHistIndex,uniqueDirs, corrDelay,PATH,FILE,bootstp_num,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the circular variance and calculate p value or significant test
% Borghuis B.G. (2003). The motion reverse correlation (MRC) method: A
% linear systems approach in the motion domain. Journal of Neuroscience
% Method 123: 153-166.
clear CV Rate Resp;
bin=0.005;
x_bin=0:bin:1;
%perm_num=1000;
figure(6);clf;orient landscape;figure(7);clf;orient landscape;figure(8);clf;orient landscape;
for (Xpatch = 1:numXpatches)
    for (Ypatch = 1:numYpatches) 
        for (DelayIndex = 1:spikeHistIndex)  
            clear Response; Response=squeeze(g_spikeTracker_f(Xpatch,Ypatch,DelayIndex, uniqueDirs+1))';
%             DirectionNumber=squeeze(g_dirCounter_f(Xpatch,Ypatch,uniqueDirs+1))';
%             Response=Response./DirectionNumber;
            Resp_f=sum(Response.*exp(2*i*uniqueDirs/180*pi))/sum(Response);
            CV_f(Xpatch,Ypatch,DelayIndex)=1-abs(Resp_f);             

            [Azi_f(Xpatch,Ypatch,DelayIndex), Ele_f(Xpatch,Ypatch,DelayIndex), Amp_f(Xpatch,Ypatch,DelayIndex)] = vectorsum(Response/max(Response));
            r_f(Xpatch,Ypatch,DelayIndex) = 1-abs(Amp_f(Xpatch,Ypatch,DelayIndex)) / sum(Response/max(Response));
        end

        figure(6);subplot(numXpatches, numYpatches, Xpatch + (Ypatch-1)*numXpatches);hist(CV_f(Xpatch,Ypatch,:),x_bin);axis([0.8 1 0 100]);title('Cirular Variance');
        figure(7);subplot(numXpatches, numYpatches, Xpatch + (Ypatch-1)*numXpatches);hist(r_f(Xpatch,Ypatch,:),x_bin);title('r Variance');axis([0.6 1 0 50]);
        
        MeanCV(Xpatch,Ypatch)=mean(CV_f(Xpatch,Ypatch,:));%axis([0.8 1 0 10]);text(0.9,8,['mean=', num2str(mean(CV))]);        
%         hist_cv=hist(CV_f(Xpatch,Ypatch,:),x_bin);
        hist_r=hist(r_f(Xpatch,Ypatch,:),x_bin);
        hist_Vectorsum=hist(Amp_f(Xpatch,Ypatch,:),x_bin);
        
%         hist_cv=hist(CV(Xpatch,Ypatch,:),x_bin);
        bin_sum=0;
        n=0;
%         while(n<=(CV_Max(Xpatch,Ypatch)/bin))
%             n=n+1;
%             bin_sum=bin_sum+hist_cv(1,n);
%             %p(Xpatch,Ypatch)=(perm_num-bin_sum)/perm_num;
%             CVp(Xpatch,Ypatch)=bin_sum/sum(hist_cv);
%         end
        
        bin_sum=0;
        n=0;
        while (n<=(r_Max(Xpatch,Ypatch)/bin))
            n=n+1;
            bin_sum=bin_sum+hist_r(1,n);
            rp(Xpatch,Ypatch)=bin_sum/sum(hist_r);
        end      
        
        bin_sum=0;
        n=0;
        while (n<=(Amp_Max(Xpatch,Ypatch)/bin) & n<200)
            n=n+1;
            bin_sum=bin_sum+hist_Vectorsum(1,n);
            VectSp(Xpatch,Ypatch)=1-bin_sum/sum(hist_Vectorsum);
        end
        
        
%         bin_sum = 0;
%         n = 0; 
%         while ( bin_sum < 0.025*sum(hist_cv) )   % define confidential value to be 0.05, now consider one side only which is 0.025 of confidence
%             n = n+1;
%             bin_sum = bin_sum + hist_cv(1, n);      
%             pp(Xpatch,Ypatch) = CV_Max(Xpatch,Ypatch) - n * bin ;    % calculate what HTI value is thought to be significant different            
%         end
        figure(4);subplot(numXpatches,numYpatches, Xpatch + (Ypatch-1)*numXpatches); 
        x_lim=xlim;y_lim=ylim;
%         text(0.5*x_lim(2),0.9*y_lim(2),['CirVar=',num2str(CV_Max(Xpatch,Ypatch))],'FontSize', 8);
%         text(0.5*x_lim(2),0.75*y_lim(2),['p=',num2str(VectSp(Xpatch,Ypatch))],'FontSize', 8);        
        text(0.5*x_lim(2),0.90*y_lim(2),['VectAmp=',num2str(Amp_Max(Xpatch,Ypatch))],'FontSize', 8);        
%         text(0.5*x_lim(2),0.75*y_lim(2),['r=',num2str(r_Max(Xpatch,Ypatch))],'FontSize', 8);        
        text(0.5*x_lim(2),0.70*y_lim(2),['rp=',num2str(rp(Xpatch,Ypatch))],'FontSize', 8);  
        text(0.5*x_lim(2),0.50*y_lim(2),['VectSump=',num2str(VectSp(Xpatch,Ypatch))],'FontSize', 8); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the Spatial_RF map with a two-dimensional Gaussian function
step_x=FieldSize_x/(numXpatches-1);
step_y=FieldSize_y/(numYpatches-1);
x_interp0=FieldLoc_x-0.5*FieldSize_x:step_x:FieldLoc_x+0.5*FieldSize_x;
y_interp0=FieldLoc_y-0.5*FieldSize_y:step_y:FieldLoc_y+0.5*FieldSize_y;
x_interp=[fliplr(FieldLoc_x-0.5*FieldSize_x:-step_x:floor(-1*(View_Dist+step_x)/step_x)*step_x) x_interp0(2:end-1) FieldLoc_x+0.5*FieldSize_x:step_x:ceil((View_Dist+step_x)/step_x)*step_x];
y_interp=[fliplr(FieldLoc_y-0.5*FieldSize_y:-step_y:floor(-1*(View_Dist+step_y)/step_y)*step_y) y_interp0(2:end-1) FieldLoc_y+0.5*FieldSize_y:step_y:ceil((View_Dist+step_x)/step_y)*step_y];
% figure;clf;orient landscape;contourf(x_interp0,y_interp0,flipud(cplot_wrapped));axis([-1*View_Dist View_Dist -1*View_Dist View_Dist]);
minLength=min([length(x_interp) length(y_interp)]);
x_interp=x_interp(1:minLength);
y_interp=y_interp(1:minLength);
Start_X=FieldLoc_x-0.5*FieldSize_x;StartpatchX=find(x_interp==Start_X);%x_interp(StartpatchX)
Start_Y=FieldLoc_y+0.5*FieldSize_y;StartpatchY=find(fliplr(y_interp)==Start_Y);
%%%%%%%%%%%
%plot RF from raw data
z_gauss1 = zeros(length(x_interp), length(y_interp))+min(min(cplot_wrapped));
for Xpatch=1:numXpatches
    for Ypatch=1:numYpatches       
        z_gauss1(StartpatchY+Ypatch-1,StartpatchX+Xpatch-1)=cplot_wrapped(Ypatch,Xpatch);        
    end
end
figure(8);clf;subplot(2,2,1);contourf(x_interp,y_interp, flipud(z_gauss1));colorbar;grid;axis([-1*View_Dist View_Dist -1*View_Dist View_Dist]);
clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'XTick',x);
ticks = atan(x/View_Dist)/pi*180;
set(gca,'XTickLabel', ticks);
clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'YTick',y);
ticks = atan(y/View_Dist)/pi*180;
set(gca,'yTickLabel', ticks);
clear y_lim y_range;
y_lim = ylim;
y_range = y_lim(2)-y_lim(1);
xlabel('X Direction (deg)');ylabel('Y Direction (deg)');Title=['Time= -',num2str(maxDelay),'ms'];title(Title);
text(x_lim(1), y_lim(2)+.15*y_range, [PATH, FILE]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate the goodness of the spatial RF
tempmatrix=cplot_wrapped;
Distance_Ori=DistanceEvaluate1(tempmatrix,VectSp);
%do the permutation
for k=1:1000
    ind=randperm(size(tempmatrix,1)*size(tempmatrix,2));
    matrixind=reshape(ind,size(tempmatrix,1),size(tempmatrix,2));
    permutematrix=tempmatrix(matrixind);
    Distance_Permute(k)=DistanceEvaluate(permutematrix);
end
%test the significance
bin=0.005;
x_bin=0:bin:1;
hist_Dist=hist(Distance_Permute/max(Distance_Permute),x_bin);
bin_sum=0;
n=0;
while(n<=((Distance_Ori/max(Distance_Permute))/bin))
    n=n+1;
    bin_sum=bin_sum+hist_Dist(1,n);   
    Dist_p=bin_sum/sum(hist_Dist);
end
%%%%%%%%%%%%%%
%Fitted the RF with 2D Gaussian
z_gauss0=flipud(cplot_wrapped);
for Xpatch=1:length(x_interp0)
    for Ypatch=1:length(y_interp0)
        temp(((Xpatch-1)*(length(y_interp0)) + Ypatch), 1) = x_interp0(Xpatch);
        temp(((Xpatch-1)*(length(y_interp0)) + Ypatch), 2) = y_interp0(Ypatch);
        temp(((Xpatch-1)*(length(y_interp0)) + Ypatch), 3) = z_gauss0(Ypatch, Xpatch);
    end
end
% x_interp1=x_interp0(1)+FieldSize_x/(2*numXpatches): FieldSize_x/numXpatches:x_interp0(1,end);
% y_interp1=y_interp0(1)+FieldSize_y/(2*numYpatches): FieldSize_y/numYpatches:y_interp0(1,end);
% for Xpatch=1:length(x_interp1)
%     for Ypatch=1:length(y_interp1)
%         temp(((Xpatch-1)*(length(y_interp1)) + Ypatch), 1) = x_interp1(Xpatch);
%         temp(((Xpatch-1)*(length(y_interp1)) + Ypatch), 2) = y_interp1(Ypatch);
%         temp(((Xpatch-1)*(length(y_interp1)) + Ypatch), 3) = z_gauss0(Ypatch, Xpatch);
%     end
% end
means(:,1) = temp(:,1);
means(:,2) = temp(:,2);
means(:,3) = temp(:,3);
raw=means;
mean_graph = zeros(length(x_interp0), length(y_interp0));
for Xpatch=1:size(mean_graph,1)
    mean_graph(:,Xpatch)=means((Xpatch-1)*size(mean_graph,2)+1:Xpatch*size(mean_graph,2),3);
end

pars = gauss2Dfit(means,raw);%2D_Gaussion Fit
for Xpatch=1:length(x_interp0)
    for Ypatch = 1:length(y_interp0)        
        Matrix_2DGaussian(Xpatch,Ypatch) =  gauss2Dfunc(x_interp0(Xpatch),y_interp0(Ypatch), pars);
    end
end
Matrix_2DGaussian=rot90(Matrix_2DGaussian);
z_gauss2 = zeros(length(x_interp), length(y_interp))+min(min(cplot_wrapped));
for Xpatch=1:numXpatches
    for Ypatch=1:numYpatches       
        z_gauss2(StartpatchY+Ypatch-1,StartpatchX+Xpatch-1)=Matrix_2DGaussian(Ypatch,Xpatch);        
    end
end
[r,p2]=corrcoef(cplot_wrapped,Matrix_2DGaussian);%[r,p2]=corrcoef(mean_graph,z_gauss);
R_Value=r(1,2);
P_Value=p2(1,2);clear r p2;
figure(8);subplot(2,2,2);contourf(x_interp,y_interp,flipud(z_gauss2));colorbar;grid;axis([-1*View_Dist View_Dist -1*View_Dist View_Dist]);
clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'XTick',x);
ticks = atan(x/View_Dist)/pi*180;
set(gca,'XTickLabel', ticks);
clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'YTick',y);
ticks = atan(y/View_Dist)/pi*180;
set(gca,'yTickLabel', ticks);
xlabel('X Direction (deg)');ylabel('Y Direction (deg)');title('2D Gaussian Fitted RF');
%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation
x_interp_new = -1*View_Dist: 0.2 : View_Dist;
y_interp_new = -1*View_Dist: 0.2 : View_Dist;

for Xpatch=1:length(x_interp_new)
    for Ypatch = 1:length(y_interp_new)
        z_gauss_Interpolation(Xpatch,Ypatch) =  gauss2Dfunc(x_interp_new(Xpatch),y_interp_new(Ypatch), pars);
    end
end
zz=flipud(rot90(z_gauss_Interpolation));
%figure(8);subplot(2,2,3);contourf(x_interp_new,y_interp_new,z_gauss_new);colorbar;grid
figure(8);subplot(2,2,3);contourf(x_interp_new,y_interp_new,zz);colorbar;grid
clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'XTick',x);
ticks = atan(x/View_Dist)/pi*180;
set(gca,'XTickLabel', ticks);
clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
set(gca,'YTick',y);
ticks = atan(y/View_Dist)/pi*180;
set(gca,'yTickLabel', ticks);
xlabel('X Direction (deg)');ylabel('Y Direction (deg)');title('RF with interpolation'); 

x_lim =xlim;y_lim = ylim;
string=sprintf('Base Rate =%1.3f', pars(1));
text(1.8*x_lim(1),y_lim(2)+10,string,'FontSize',8);
string=sprintf('Amplitude =%1.3f', pars(2));
text(1.8*x_lim(1),y_lim(2)+5,string,'FontSize',8);
string=sprintf('FWHM_X =%1.3f', atan(2.35*pars(4)/View_Dist)/pi*180);%(FWHM (full width half maximum) by the equation: sigma=FWHM/(2*SQRT(ALOG(2))), this is typically written as: sigma=FWHM/2.35; Peter Neri. Spatial Integration of Optic flow signals in fly motion-sensitive neurons. Neurophysiol. 95:1608-1619,2006)
text(1.8*x_lim(1),y_lim(2),string,'FontSize',8);
string=sprintf('FWHM_Y =%1.3f', atan(2.35*pars(6)/View_Dist)/pi*180);
text(1.8*x_lim(1),y_lim(2)-5,string,'FontSize',8);
string=sprintf('Fitted CTR =(%1.3f,%1.3f)',atan(pars(3)/View_Dist)/pi*180,atan(pars(5)/View_Dist)/pi*180);
text(1.8*x_lim(1),y_lim(2)-10,string,'FontSize',8);
string=sprintf('Eccentricity =%1.3f',atan(sqrt(pars(3).^2+ pars(5).^2)/View_Dist)/pi*180);
text(1.8*x_lim(1),y_lim(2)-15,string,'FontSize',8);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fit the Spatial_RF map with a two-dimensional Gaussian function
% step_x=FieldSize_x/(numXpatches-1);
% step_y=FieldSize_y/(numYpatches-1);
% x_interp0=FieldLoc_x-0.5*FieldSize_x:step_x:FieldLoc_x+0.5*FieldSize_x;
% y_interp0=FieldLoc_y-0.5*FieldSize_y:step_y:FieldLoc_y+0.5*FieldSize_y;
% x_interp=[fliplr(FieldLoc_x-0.5*FieldSize_x:-step_x:floor(-1*View_Dist/step_x)*step_x) x_interp0(2:end-1) FieldLoc_x+0.5*FieldSize_x:step_x:ceil(View_Dist/step_x)*step_x];
% y_interp=[fliplr(FieldLoc_y-0.5*FieldSize_y:-step_y:floor(-1*View_Dist/step_y)*step_y) y_interp0(2:end-1) FieldLoc_y+0.5*FieldSize_y:step_y:ceil(View_Dist/step_y)*step_y];
% % figure;clf;orient landscape;contourf(x_interp0,y_interp0,flipud(cplot_wrapped));axis([-1*View_Dist View_Dist -1*View_Dist View_Dist]);
% 
% minLength=min([length(x_interp) length(y_interp)]);
% x_interp=x_interp(1:minLength);
% y_interp=y_interp(1:minLength);
% z_gauss = zeros(length(x_interp), length(y_interp));
% 
% Start_X=FieldLoc_x-0.5*FieldSize_x;StartpatchX=find(x_interp==Start_X);%x_interp(StartpatchX)
% Start_Y=FieldLoc_y+0.5*FieldSize_y;StartpatchY=find(fliplr(y_interp)==Start_Y);
% 
% for Xpatch=1:numXpatches
%     for Ypatch=1:numYpatches       
%         z_gauss(StartpatchY+Ypatch-1,StartpatchX+Xpatch-1)=cplot_wrapped(Ypatch,Xpatch);        
%     end
% end
% 
% %figure;orient landscape;contourf(x_interp,y_interp, flipud(z_gauss));grid% figure(8);clf;orient landscape;contourf(x_interp,y_interp,z_gauss);
% z_gauss0=flipud(z_gauss);
% for Xpatch=1:length(x_interp)
%     for Ypatch=1:length(y_interp)
%         temp(((Xpatch-1)*(length(y_interp)) + Ypatch), 1) = x_interp(Xpatch);
%         temp(((Xpatch-1)*(length(y_interp)) + Ypatch), 2) = y_interp(Ypatch);
%         temp(((Xpatch-1)*(length(y_interp)) + Ypatch), 3) = z_gauss0(Ypatch, Xpatch); %z_gauss(Xpatch, Ypatch);     
%     end
% end
% means(:,1) = temp(:,1);
% means(:,2) = temp(:,2);
% means(:,3) = temp(:,3);
% raw=means;
% mean_graph = zeros(length(x_interp), length(y_interp));
% for Xpatch=1:size(mean_graph,1)
%     mean_graph(:,Xpatch)=means((Xpatch-1)*size(mean_graph,2)+1:Xpatch*size(mean_graph,2),3);
% end
% figure(8);clf;subplot(2,2,1);contourf(x_interp, y_interp, mean_graph);colorbar;%grid;%axis([-1*View_Dist View_Dist -1*View_Dist View_Dist]);
% clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'XTick',x);
% ticks = atan(x/View_Dist)/pi*180;
% set(gca,'XTickLabel', ticks);
% clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'YTick',y);
% ticks = atan(y/View_Dist)/pi*180;
% set(gca,'yTickLabel', ticks);
% clear y_lim y_range;
% y_lim = ylim;
% y_range = y_lim(2)-y_lim(1);
% xlabel('X Direction (deg)');ylabel('Y Direction (deg)');Title=['Time= -',num2str(maxDelay),'ms'];title(Title);
% text(x_lim(1), y_lim(2)+.15*y_range, [PATH, FILE]); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Evaluate the goodness of the spatial RF
% tempmatrix=cplot_wrapped;
% Distance_Ori=DistanceEvaluate1(tempmatrix,VectSp);
% %do the permutation
% for k=1:1000
%     ind=randperm(size(tempmatrix,1)*size(tempmatrix,2));
%     matrixind=reshape(ind,size(tempmatrix,1),size(tempmatrix,2));
%     permutematrix=tempmatrix(matrixind);
%     VectSp_matrix=VectSp(matrixind);
%     Distance_Permute(k)=DistanceEvaluate1(permutematrix,VectSp_matrix);  
% end
% %test the significance
% bin=0.005;
% x_bin=0:bin:1;
% hist_Dist=hist(Distance_Permute/max(Distance_Permute),x_bin);
% 
% bin_sum=0;
% n=0;
% while(n<=((Distance_Ori/max(Distance_Permute))/bin))
%     n=n+1;
%     bin_sum=bin_sum+hist_Dist(1,n);   
%     Dist_p=bin_sum/sum(hist_Dist);
% end
% 
% %%%%%%%%%%%%%%
% %2D_Gaussion Fit
% pars = gauss2Dfit(means,raw);
% unique_x_pos = -1*View_Dist:step_x:View_Dist;
% unique_y_pos = -1*View_Dist:step_y:View_Dist;
% x_interp_new = unique_x_pos(1): 0.2 : unique_x_pos(length(unique_x_pos));
% y_interp_new = unique_y_pos(1): 0.2 : unique_y_pos(length(unique_y_pos));
% z_gauss = zeros(length(x_interp), length(y_interp));
% 
% for Xpatch=1:length(x_interp)
%     for Ypatch = 1:length(y_interp)        
%         z_gauss(Xpatch,Ypatch) =  gauss2Dfunc(x_interp(Xpatch),y_interp(Ypatch), pars);
%     end
% end
% 
% for Xpatch=1:length(x_interp_new)
%     for Ypatch = 1:length(y_interp_new)
%         z_gauss_new(Xpatch,Ypatch) =  gauss2Dfunc(x_interp_new(Xpatch),y_interp_new(Ypatch), pars);
%     end
% end
% [r,p2]=corrcoef(mean_graph,flipud(rot90(z_gauss)));%[r,p2]=corrcoef(mean_graph,z_gauss);
% R_Value=r(1,2);
% P_Value=p2(1,2);clear r p2;
% figure(8);subplot(2,2,2);contourf(x_interp,y_interp,flipud(rot90(z_gauss)));colorbar;%grid
% clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'XTick',x);
% ticks = atan(x/View_Dist)/pi*180;
% set(gca,'XTickLabel', ticks);
% clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'YTick',y);
% ticks = atan(y/View_Dist)/pi*180;
% set(gca,'yTickLabel', ticks);
% xlabel('X Direction (deg)');ylabel('Y Direction (deg)');title('2D Gaussian Fitted RF');
% zz=flipud(rot90(z_gauss_new));
% %figure(8);subplot(2,2,3);contourf(x_interp_new,y_interp_new,z_gauss_new);colorbar;grid
% figure(8);subplot(2,2,3);contourf(x_interp_new,y_interp_new,zz);colorbar;%grid
% 
% rf_xctr=data.one_time_params(RF_XCTR);rf_yctr=data.one_time_params(RF_YCTR);
% x_lim =xlim;y_lim = ylim;
% string=sprintf('Base Rate =%1.3f', pars(1));
% text(1.8*x_lim(1),y_lim(2)+3,string,'FontSize',8);
% string=sprintf('Amplitude =%1.3f', pars(2));
% text(1.8*x_lim(1),y_lim(2)-2,string,'FontSize',8);
% string=sprintf('FWHM_X =%1.3f', atan(2.35*pars(4)/View_Dist)/pi*180);%(FWHM (full width half maximum) by the equation: sigma=FWHM/(2*SQRT(ALOG(2))), this is typically written as: sigma=FWHM/2.35; Peter Neri. Spatial Integration of Optic flow signals in fly motion-sensitive neurons. Ypatch.Neurophysiol. 95:1608-1619,2006)
% text(1.8*x_lim(1),y_lim(2)-7,string,'FontSize',8);
% string=sprintf('FWHM_Y =%1.3f', atan(2.35*pars(6)/View_Dist)/pi*180);
% text(1.8*x_lim(1),y_lim(2)-12,string,'FontSize',8);
% string=sprintf('Fitted CTR =(%1.3f,%1.3f)',atan(pars(3)/View_Dist)/pi*180,atan(pars(5)/View_Dist)/pi*180);
% text(1.8*x_lim(1),y_lim(2)-17,string,'FontSize',8);
% % string=sprintf('Eccentricity =%1.3f',atan(sqrt(pars(3).^2+ pars(5).^2)/View_Dist)/pi*180);
% % text(1.8*x_lim(1),y_lim(2)-22,string,'FontSize',8);
% 
% clear x;x=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'XTick',x);
% ticks = atan(x/View_Dist)/pi*180;
% set(gca,'XTickLabel', ticks);
% clear y;y=[-1 -0.4142 0 0.4142 1]*View_Dist;
% set(gca,'YTick',y);
% ticks = atan(y/View_Dist)/pi*180;
% set(gca,'yTickLabel', ticks);
% xlabel('X Direction (deg)');ylabel('Y Direction (deg)');title('RF with interpolation'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SpikeChan==1
    Directory=['save C:\Aihua\z_TempOutputs\SU']
else
    Directory=['save C:\Aihua\z_TempOutputs\MU']
end
save_name = 'cplot'; % name of the  mat file with all the Original CrossCorrelation Histrogram 
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= cplot_wrapped;'])
save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'Matrix_2DGaussian'; % name of the  mat file with all the Original CrossCorrelation Histrogram 
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= Matrix_2DGaussian;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'z_gauss_interp'; % name of the  mat file with all the Original CrossCorrelation Histrogram
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= z_gauss_Interpolation;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'R_CurveFit'; % name of the mat
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= R_CurveFit;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'Dist_p'; % name of the mat
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= Dist_p;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'CirVar'; % name of the mat of circular variance from the noise (Reverse Correlaiton_forward)
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= CV_f;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'CirVar_Max'; % name of the mat of circular variance from the noise (Reverse Correlaiton_forward)
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= CV_Max;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'VectSum_p'; % 
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= VectSp;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'VectAmp'; % name of the mat of circular variance from the noise (Reverse Correlaiton_forward)
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= Amp_Max;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'r_Max'; % name of the mat of circular variance from the noise (Reverse Correlaiton_forward)
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= r_Max;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'r_f'; % name of the mat of circular variance from the noise (Reverse Correlaiton_forward)
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= r_f;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

save_name = 'r_p'; % tuning index
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= rp;'])
clear save_file;save_file=[Directory(6:end) save_name];
if (exist([save_file '.mat'], 'file') == 0)    %file does not yet exist  
    eval([Directory save_name spacer dum_name]);
else
    eval([Directory save_name spacer dum_name ' -APPEND']);   
end

%--------------------------------------------------------------------------
%Write out some summary data to a cumulative summary file
sprint_txt=['%s\t %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %1.3f\t   %6.3f\t   %6.3f\t    %6.3f\t   %6.3f\t   %6.3f\t'];
if SpikeChan==1
    FILEName=[FILE '_SU'];
else
    FILEName=[FILE '_MU'];
end

%buff=sprintf(sprint_txt, FILEName, pars(1), pars(2), atan(2.35*pars(4)/View_Dist)/pi*180, atan(2.35*pars(6)/View_Dist)/pi*180, atan(pars(3)/View_Dist)/pi*180,atan(pars(5)/View_Dist)/pi*180, rf_xctr, rf_yctr, R_Value, P_Value, numXpatches, numYpatches,Dist_p);
buff=sprintf(sprint_txt, FILEName, pars(1), pars(2), pars(3),pars(4),pars(5),pars(6), View_Dist,atan(2.35*pars(4)/View_Dist)/pi*180, atan(2.35*pars(6)/View_Dist)/pi*180, atan(pars(3)/View_Dist)/pi*180,atan(pars(5)/View_Dist)/pi*180,atan(sqrt(pars(3).^2+pars(5).^2)/View_Dist)/pi*180,R_Value, P_Value, numXpatches, numYpatches,Dist_p);
outfile =['C:\Aihua\z_TempOutputs\ReverseCorr.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Base Rate\t Amplitude\t XCenter_cm\t SigmaX_cm\t YCenter_cm\t SigmaY_cm\t ViewDistance\t FWHMx\t FWHMy\t Fitted_CTRx\t Fitted_CTRy\t Eccentricity\t R_Value\t P_Value\t numXpatches\t numYpatches\t Dist_p\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);

%---------------------------------------------------------------------------------------

return;