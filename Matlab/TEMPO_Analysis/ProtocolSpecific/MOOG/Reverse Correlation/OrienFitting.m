function [F, x, x0] =OrienFitting(Input)
%[F, x, x0] =OrienFitting(Input)
%Responses to different directions were fit with a wrapped Gaussian
%Equation: 
global xData yData;
yData=Input; clear Input
xData=[0    45    90   135   180   225   270   315 ]*pi/180;
lb=[0.001 -2*pi pi/6 0];
ub=[1.5*(max(yData)-min(yData)) 2*pi 2*pi 0.8*max(yData)];
param_label=['A          mu             sigma           DC'];%A: tuning amplitude // mu: preferred direction // sigma: bandwidth //DC: baseline firing rate
x0(1)=max(yData)-min(yData);
find(yData==max(yData));
x0(2)=xData(find(yData==max(yData)));
x0(4)=min(yData);
%search for best starting values of x0(3)
N=30;
min_err=9999999999999999.99;%???????
x3range = lb(3) : (ub(3)-lb(3))/N : ub(3);
for i=1:N 
    x3temp=x3range(i);
    if x3temp==0
        x3temp=0.0001;
    end
    x_temp=[x0(1) x0(2) x3temp x0(4)];
    error = DirCurvefit_err(x_temp);
%          yfit=DirCurvefit(x_temp,xData);
%          error = norm(sqrt(yfit)-sqrt(yData))^2;
    if (error<min_err)
        x3min=x3temp;
        min_err=error;
    end    
end
x0(3)=x3min;
F=DirCurvefit(x0,xData);
% figure;plot(xData,yData,'bo');
% hold on;plot(xData,F,'r');

options = optimset('MaxFunEvals', 10000, 'MaxIter', 5000, 'LargeScale', 'off', 'LevenbergMarquardt', 'on', 'Display', 'off');
A = []; b = []; Aeq = []; beq = []; nonlcon = [];

% fit multiple times with some jitter in the initial params
N_reps = 30;
wiggle = 0.3;
clear testpars fval_temp exitflag_temp;
for j=1:N_reps  
    rand_factor = rand(length(x0),1) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_x0 = x0' .* rand_factor;
    [testpars{j}, fval_temp(j), exitflag_temp(j)] = fmincon('DirCurvefit_err', temp_x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    err(j) = DirCurvefit_err(testpars{j});
    yfit_test{j} = (DirCurvefit(testpars{j},xData));
    [coef_test,P_test] = corrcoef(yData',yfit_test{j}');   % R^2's using means
    rsquared_test(j) = coef_test(1,2)^2;
end

% so the best-fit values of x are:
best_err = find(err == min(err));
best_r2 = find(rsquared_test(best_err) == max(rsquared_test(best_err)));  % if multiple occurences of lowest 
if length(best_r2)>1
    best_r3=best_r2(1);clear best_r2;
    best_r2=best_r3;clear best_r3;
end
x = testpars{best_err(best_r2)};                                     % error, use the one with best R^2
% and the rest...
final_error = min(err);
fval = fval_temp(best_err(best_r2));
exitflag = exitflag_temp(best_err(best_r2));
yfit = (DirCurvefit(x,xData));
[coef,P] = corrcoef(yData,yfit);%[coef,P] = corrcoef(mean(yData)',mean(yfit)');   % R^2's using means
rsquared = coef(1,2)^2;
p_fit = P(1,2);
x_0 = x0;

%finish plotting
% x_smooth = 0:0.01:2*pi;
% y_smooth = (DirCurvefit(x,x_smooth));
% xdata_tran = [0    45    90   135   180   225   270   315]*pi/180;%xdata_tran = [-90 -45 0 45 90 135 180 225 270] * pi/180;
% x_smooth_tran = xdata_tran(1):0.01:xdata_tran(end);
% y_smooth_tran = (DirCurvefit(x,x_smooth_tran));
% plot(x_smooth_tran,y_smooth_tran,'b');
% xlim( [xdata_tran(1) xdata_tran(end)] );
% set(gca, 'xtick', xdata_tran);
% %set(gca, 'xdir' , 'reverse');
% set(gca, 'xticklabel','0|45|90|135|180|225|270|315');
% xlabel('Azimuth');
% legend('Original Data','Initial Fitting' ,'Jitter Fitting',2)

% % show params for each fit
% param_text = [num2str(x(1),4) '    ' num2str(x(2)*180/pi,4) '    ' num2str(x(3)*180/pi,4) '    ' num2str(x(4),4)   ];  
% x0_text = [num2str(x0(1),4) '    ' num2str(x0(2)*180/pi,4) '    ' num2str(x0(3)*180/pi,4) '    ' num2str(x0(4),4) ];
% y_lim = ylim;
% y_range = y_lim(2)-y_lim(1);
% text(4, y_lim(2)+.20*y_range, param_label);
% text(4, y_lim(2)+.12*y_range, x0_text);
% text(4, y_lim(2)+.04*y_range, param_text);
% 
% % take peak of fitted function as preferred heading 
% peak = x_smooth(find(y_smooth == max(y_smooth))) * 180/pi;