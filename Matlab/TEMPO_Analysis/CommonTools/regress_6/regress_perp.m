function result =regress_perp(x,y,option,alpha)
if isempty([x(:);y(:)])
    result = [];
    return
end

% regress_perp.m 
% 
% [slope, intercept]=regress(x,y);
% [slope, intercept]=regress(x,y,alpha);
% [slope, intercept, slopeInt, interceptInt, r, p]=regress_perp(x,y,option)
% [slope, intercept, slopeInt, interceptInt, r, p]=regress_perp(x,y,option,alpha)
%
% x and y are two column vectors that contain the data set's x and y value, the
% two's size must match. The function return slope, intercept, slopes 95%
% confidence intervel, intercept 95% confidence intervel, r value,
% and p value, the user may specify the confidence interval by input alpha.
% (alpha input is optional)
% option:
% 1--free intercept (default)
% 2--set intercept 0
% 3--set slope 2 (???)
% 
% this script estimates the linear regression line by minimizing the 
% perpendicular offset. This is different from the general linear regression
% method which minimizes the vertical offset. (See function regress.m)
% disp('***** b--slope; a--intercept; bint--slope confidence int; aint--inter confidence int *****')
% disp('***** 1--free intercept (default); 2--set intercept 0*****')
if  nargin < 2,              
    error('REGRESS_PERP requires at least two input arguments.');      
elseif nargin == 2
    option = 1;
    alpha = 0.05;
elseif nargin == 3
    alpha = 0.05;
end

if option~=1 & option~=2 & option~=3
    error('option input is not recognized');
end

meanX = mean(x); meanY = mean(y); 

%estimating slope and intercept via nonlinear least square fitting
%the fitting minimizes perpendicular offset (point2line distance)
if option==1 
    
    %     result.k = regress_slope([x y]);
    %     result.b = regress_int([x y]);
    
    % HH20180613
    coeff = pca([x y]);
    result.k = coeff(2) / coeff(1);
    result.b = meanY- result.k * meanX;
    
elseif option==2  % fixed offset
    
    result.k = regress_slope_int([x y]); 
    result.b = 0;
    
else  % fixed slope
    result.k = 2;
    result.b = regress_int_slope([x y]);
    
end
    
if isnan(result.k) | isnan(result.b)
    error('The program failed to extact slope and intercept from the data set via least square');
end

opt = statset('UseParallel',false);

if nargin > 3
    % clear index;
    %estimating confidence interval for slope
    if option==1
        
        % HH20180613
        [bootPCA, bootSample] = bootstrp(1000, @pca, [x y],'Options',opt);
        bootMeanX = mean(x(bootSample),1)';
        bootMeanY = mean(y(bootSample),1)';
        
        result.bootK = bootPCA(:,2) ./ bootPCA(:,1);
        result.bootB = bootMeanY - result.bootK .* bootMeanX;
        
        % slope_sam=bootstrp(1000,'regress_slope',[x y],'Options',opt);           %re-sample data via bootstrap to estimate distribution for slope

    elseif option==2
        result.bootK=bootstrp(1000,'regress_slope_int',[x y],'Options',opt);           %re-sample data via bootstrap to estimate distribution for slope
    else
        result.bootK(1:1000) = 1;
        result.bootB = bootstrp(1000,'regress_int_slope',[x y]); 
    end
    
    result.bootK=sort(result.bootK);
    result.bootB=sort(result.bootB);

    %     index=find(isnan(result.bootK));
    %
    %     if ~isempty(index)
    %         result.bootK=result.bootK(1:index(1));
    %     end
    
    result.kInterval = [result.bootK(round(alpha/2*length(result.bootK))), result.bootK(round((1-alpha/2)*length(result.bootK)))];
    result.bInterval = [result.bootB(round(alpha/2*length(result.bootB))), result.bootB(round((1-alpha/2)*length(result.bootB)))];

    % clear index;
    
    %     if option==3
    %         %         result.bootB=bootstrp(1000,'regress_int',[x y]);               % Stupid to do this twice... HH
    %         %     elseif option==2
    %         %         result.bootB=0;
    %         %     else
    %
    %     end
    
    %     index=find(isnan(result.bootB));
    %     if ~isempty(index)
    %         result.bootB=result.bootB(1:index(1));   % What's the logic???
    %     end
    %
    %     clear index;
    
    [r,p]=corrcoef(x,y);          
    result.r = r(2,1);
    result.p = p(2,1);
end

return