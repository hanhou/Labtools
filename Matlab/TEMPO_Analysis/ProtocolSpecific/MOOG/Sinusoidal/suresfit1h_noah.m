function [out]=suresfit1h_noah(frqf,x,y);

%SURESFIT1h.m 
%
%  fitting interface for single unit (su)/eyemov routines
%
% function that fits a single sinusoid to data
% to average response and stimulus cycles
% outputs: gain, phase, dc
% estimation of gain (expressed relative to ampl=peak stimulus velocity)
% phase (in deg) and dc (in deg/s) for each cycle

clear out1 out2 out3 out4 out5 tmp c

[out1,resnorm,residual1,exitflag1,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 0 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out2,resnorm,residual2,exitflag2,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 45 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out3,resnorm,residual3,exitflag3,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 90 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out4,resnorm,residual4,exitflag4,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 135 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out5,resnorm,residual5,exitflag5,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 180 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting

if exitflag1==0
    residual1=NaN;
elseif exitflag2==0
    residual2=NaN;
elseif exitflag3==0
    residual3=NaN;
elseif exitflag4==0
    residual4=NaN;
elseif exitflag5==0
    residual5=NaN;
end
tmp=[sum(residual1),sum(residual2),sum(residual3),sum(residual4),sum(residual5)];
c=find(tmp==min(tmp));
if ~isempty(c)
    c=c(1);
    eval(['out=out',num2str(c),';']);
else
    out=[NaN NaN NaN];
end

clear out1 out2 out3 out4 out5 tmp c
if isnan(out)
    out=[0 0 0];
end
[out1,resnorm,residual1,exitflag1,output,lambda,jacobian] = lsqnonlin(@sine_noah,[0 out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out2,resnorm,residual2,exitflag2,output,lambda,jacobian] = lsqnonlin(@sine_noah,[50 out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out3,resnorm,residual3,exitflag3,output,lambda,jacobian] = lsqnonlin(@sine_noah,[100 out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out4,resnorm,residual4,exitflag4,output,lambda,jacobian] = lsqnonlin(@sine_noah,[150 out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
[out5,resnorm,residual5,exitflag5,output,lambda,jacobian] = lsqnonlin(@sine_noah,[200 out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting

if exitflag1==0
    residual1=NaN;
elseif exitflag2==0
    residual2=NaN;
elseif exitflag3==0
    residual3=NaN;
elseif exitflag4==0
    residual4=NaN;
elseif exitflag5==0
    residual5=NaN;
end
tmp=[sum(residual1),sum(residual2),sum(residual3),sum(residual4),sum(residual5)];
c=find(tmp==min(tmp));
if ~isempty(c)
    c=c(1);
    eval(['out=out',num2str(c),';']);
else
    out=[NaN NaN NaN];
end


% clear out1 out2 out3 out4 out5 tmp c
% if isnan(out)
%     out=[0 0 0];
% end
% [out1,resnorm,residual1,exitflag1,output,lambda,jacobian] = lsqnonlin(@sine_noah,[out(1) out(2) 0],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
% [out2,resnorm,residual2,exitflag2,output,lambda,jacobian] = lsqnonlin(@sine_noah,[out(1) out(2) 50],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
% [out3,resnorm,residual3,exitflag3,output,lambda,jacobian] = lsqnonlin(@sine_noah,[out(1) out(2) 100],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
% [out4,resnorm,residual4,exitflag4,output,lambda,jacobian] = lsqnonlin(@sine_noah,[out(1) out(2) 150],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
% [out5,resnorm,residual5,exitflag5,output,lambda,jacobian] = lsqnonlin(@sine_noah,[out(1) out(2) 200],[],[],optimset('Display','off','MaxIter',1e6),x,y,frqf);            %Least square fitting
% 
% if exitflag1==0
%     residual1=NaN;
% elseif exitflag2==0
%     residual2=NaN;
% elseif exitflag3==0
%     residual3=NaN;
% elseif exitflag4==0
%     residual4=NaN;
% elseif exitflag5==0
%     residual5=NaN;
% end
% tmp=[sum(residual1),sum(residual2),sum(residual3),sum(residual4),sum(residual5)];
% c=find(tmp==min(tmp));
% if ~isempty(c)
%     c=c(1);
%     eval(['out=out',num2str(c),';']);
% else
%     out=[NaN NaN NaN];
% end


if ~isnan(out)
    if out(1)<0
        out(1)=-out(1);
        out(2)=out(2)+180;
    end
    out(2)=qphase(out(2));		% unwraps to [-180,180] deg
    fprintf('done!');            
else
    fprintf('failed!')
end


return

