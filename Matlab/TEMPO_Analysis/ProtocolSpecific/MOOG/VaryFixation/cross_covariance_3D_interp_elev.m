% CROSS_COVARIANCE_3D_elev.m -- same as cross covariance 3D, but 
% interpolates tuning surface in elevation also.  Gives same result, 
% because shift is only horizontal.

function DI = cross_covariance_3D_elev(resp_mat, unique_stim_type, unique_condition_num, bin, method, showplots);

shift = 180;  %(doesn't work w/ anything else right now)
x = [0 45 90 135 180 225 270 315 360];
y = [-90 -45 0 45 90];
xi = 0 : bin : 359;
xi_shift = -shift : bin : 360+shift-1;

for n = 1:length(unique_stim_type)
    for k = 1:length(unique_condition_num)
        
        resp_mat_360{k,n} = squeeze(resp_mat(k+3*(n-1),:,:));  % kluge, to make the resp_mat more intuitive
        resp_mat_360{k,n}(:,end+1) = resp_mat_360{k,n}(:,1);    
        resp_mat_360{k,n}(1,:) = resp_mat_360{k,n}(1,1);
        resp_mat_360{k,n}(end,:) = resp_mat_360{k,n}(end,1);
               
        for j = 1:length(y)    
            clear temp_y;
            temp_y = resp_mat_360{k,n}(j,:);
            yi{k,n}(j,:) = interp1(x, temp_y, xi, method); % interpolate tuning surface

            % initialize expanded versions of each tuning curve, one displaced by 'shift' degrees as a starting point
            yi_expand{k,n}(j,:) = [yi{k,n}(j, length(xi)-shift+1 : length(xi))' ; yi{k,n}(j,:)' ; yi{k,n}(j, 1 : shift)'];
            yi_shift{k,n}(j,:) = [yi{k,n}(j,:)' ; yi{k,n}(j, 1 : 2*shift)'];
        end
        
        for i = 1:length(xi)  
            clear temp_x;
            temp_x = yi{k,n}(:,i);
            yi2{k,n}(:,i) = interp1(y, temp_x, [-90:bin:89] , method);
        end
        for h = 1:180
            yi2_expand{k,n}(h,:) = [yi2{k,n}(h, length(xi)-shift+1 : length(xi))' ; yi2{k,n}(h,:)' ; yi2{k,n}(h, 1 : shift)'];
            yi2_shift{k,n}(h,:) = [yi2{k,n}(h,:)' ; yi2{k,n}(h, 1 : 2*shift)'];
        end
         
    end
    
    yi_expand = yi2_expand; yi_shift = yi2_shift;
    
    yi_shift_init = yi_shift;
    
    % now find covariance between each tuning curve and the shifted version of the other two:

% %     % Cov with itself as a test:
% %     yi_shift = yi_shift_init;  % re-initialize
% %     clear cov_self;
% %     for q = 1:2*shift
% %         clear cov_temp;
% %         cov_temp = cov(yi_expand{2,n},yi_shift{2,n});
% %         cov_self(q) = cov_temp(1,2);
% %         
% %         if showplots
% %             figure(100*n);
% %             plot(xi_shift, yi_expand{2,n}, xi_shift, yi_shift{2,n});
% %             xlim([xi_shift(1) xi_shift(end)]);
% %             title([num2str(cov_self(q)) '    ' num2str(q*bin-bin-shift) '    ' num2str((q*bin-bin-shift)/unique_condition_num(3))]);
% %             pause;
% %         end
% % 
% %         % shift yi_shift one slot over
% %         for j = 1:length(y)
% %             yi_temp = yi_shift;
% %             for r = 2:length(yi_shift{2,n})
% %                 yi_shift{2,n}(j,r) = yi_temp{2,n}(j,r-1);
% %             end
% %             yi_shift{2,n}(j,1) = yi_temp{2,n}(j,end);
% %         end
% %     end
% %     cov_self = cov_self(shift/2 + 1 : length(xi) - shift/2);
% %     DI_index = find(cov_self == max(cov_self));
% %     DI_self(n) = (DI_index - shift/2 - 1) * bin / unique_condition_num(3)
    
    % Minus with Zero:
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_MZ;
    for q = 1:2*shift
        clear cov_temp;
        cov_temp = cov(yi_expand{1,n},yi_shift{2,n});
        cov_MZ(q) = cov_temp(1,2);
        
        if showplots
            figure(100*n+1);
            plot(xi_shift, yi_expand{1,n}, xi_shift, yi_shift{2,n});
            xlim([xi_shift(1) xi_shift(end)]);
            title([num2str(cov_MZ(q)) '    ' num2str(q*bin-bin-shift) '    ' num2str((q*bin-bin-shift)/unique_condition_num(3))]);
            pause;
        end
        
        % shift yi_shift one slot over
        for j = 1:180
            yi_temp = yi_shift;
            for r = 2:length(yi_shift{2,n})
                yi_shift{2,n}(j,r) = yi_temp{2,n}(j,r-1);
            end
            yi_shift{2,n}(j,1) = yi_temp{2,n}(j,end);
        end
    end
    cov_MZ = cov_MZ(shift/2 + 1 : length(xi) - shift/2);
    DI_index = find(cov_MZ == max(cov_MZ));
    DI(1,n) = (DI_index - shift/2 - 1) * bin / unique_condition_num(3);

    % Minus with Plus
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_MP;
    for q = 1:2*shift
        clear cov_temp;
        cov_temp = cov(yi_expand{1,n},yi_shift{3,n});
        cov_MP(q) = cov_temp(1,2);

        if showplots        
            figure(100*n+2);
            plot(xi_shift, yi_expand{1,n}, xi_shift, yi_shift{3,n});
            xlim([xi_shift(1) xi_shift(end)]);
            title([num2str(cov_MP(q)) '    ' num2str(q*bin-bin-shift) '    ' num2str((q*bin-bin-shift)/(unique_condition_num(3)*2))]);
            pause;
        end
        
        % shift yi_shift one slot over
        for j = 1:180
            yi_temp = yi_shift;
            for r = 2:length(yi_shift{3,n})
                yi_shift{3,n}(j,r) = yi_temp{3,n}(j,r-1);
            end
            yi_shift{3,n}(j,1) = yi_temp{3,n}(j,end);
        end
    end
    cov_MP = cov_MP(shift/2 + 1 : length(xi) - shift/2);
    DI_index = find(cov_MP == max(cov_MP));
    DI(2,n) = (DI_index - shift/2 - 1) * bin / (unique_condition_num(3)*2);
           
    % Zero with Plus
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_ZP;
    for q = 1:2*shift
        clear cov_temp;
        cov_temp = cov(yi_expand{2,n},yi_shift{3,n});
        cov_ZP(q) = cov_temp(1,2);

        if showplots
            figure(100*n+3);
            plot(xi_shift, yi_expand{2,n}, xi_shift, yi_shift{3,n});
            xlim([xi_shift(1) xi_shift(end)]);
            title([num2str(cov_ZP(q)) '    ' num2str(q*bin-bin-shift) '    ' num2str((q*bin-bin-shift)/unique_condition_num(3))]);
            pause;
        end
        
        % shift yi_shift one slot over
        for j = 1:180
            yi_temp = yi_shift;
            for r = 2:length(yi_shift{3,n})
                yi_shift{3,n}(j,r) = yi_temp{3,n}(j,r-1);
            end
            yi_shift{3,n}(j,1) = yi_temp{3,n}(j,end);
        end
    end
    cov_ZP = cov_ZP(shift/2 + 1 : length(xi) - shift/2);
    DI_index = find(cov_ZP == max(cov_ZP));
    DI(3,n) = (DI_index - shift/2 - 1) * bin / unique_condition_num(3);

end

return;