% CROSS_COVARIANCE.m -- compute covariance between heading tuning curves in a
% 1D Vary_Fixation experiment.  Curves obtained at different eye positions
% are interpolated, then shifted systematically with respect to each other.
% The shift that produces the largest covariance (normalized to the shift in 
% eye position) is taken as the shift ratio or 'Displacement Index' (DI; 
% see Avillac, Deneve, Olivier, Pouget, and Duhamel, Nat Nsci. 2005).

function DI = cross_covariance(resp_mat, unique_stim_type, unique_condition_num, bin, method, showplots, x);

shift = 180;  %(doesn't work w/ anything else right now)
xi = 0 : bin : 359;
xi_shift = -shift + x(1) : bin : x(end) + shift - 1;

resp_mat_360 = resp_mat;
for n = 1:length(unique_stim_type)

    for k = 1:length(unique_condition_num)
        resp_mat_360(k,n,length(x)) = resp_mat(k,n,1);
        clear temp_y;
        for i = 1:length(x)
            temp_y(i) = resp_mat_360(k,n,i);
        end
        
        yi{k,n} = interp1(x, temp_y, xi, method); % interpolate tuning curve

        % initialize expanded (wrapped) versions of each tuning curve, one displaced by 'shift' degrees as a starting point
        yi_expand{k,n} = [yi{k,n}(length(xi) - (shift/bin)+1 : length(xi))' ; yi{k,n}' ; yi{k,n}(1 : shift/bin)'];
        yi_shift{k,n} = [yi{k,n}' ; yi{k,n}(1 : 2*shift/bin)'];
        
    end
    
    yi_shift_init = yi_shift;
    
    % now find covariance between each tuning curve and the shifted version of the other two:

% %     % Cov a tuning curve with itself, as a sanity check:
% %     yi_shift = yi_shift_init;  % re-initialize
% %     clear cov_self;
% %     for q = 1:2*(shift/bin)
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
% %         yi_temp = yi_shift;
% %         for r = 2:length(yi_shift{2,n})
% %             yi_shift{2,n}(r) = yi_temp{2,n}(r-1);
% %         end
% %         yi_shift{2,n}(1) = yi_temp{2,n}(end);
% %     end
% %     cov_self = cov_self((shift/bin)/2 + 1 : length(xi) - (shift/bin)/2);  % limit possible DI's to only plausible values (+/- 90 deg)
% %     DI_index = find(cov_self == max(cov_self));
% %     DI_self(n) = (DI_index - (shift/bin)/2 - 1) * bin / unique_condition_num(3)
    
    % Minus with Zero:
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_MZ;
    for q = 1:2*(shift/bin)
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
        yi_temp = yi_shift;
        for r = 2:length(yi_shift{2,n})
            yi_shift{2,n}(r) = yi_temp{2,n}(r-1);
        end
        yi_shift{2,n}(1) = yi_temp{2,n}(end);
    end
    cov_MZ = cov_MZ((shift/bin)/2 + 1 : length(xi) - (shift/bin)/2); % limit possible DI's to only plausible values (+/- 90 deg)
    DI_index = find(cov_MZ == max(cov_MZ));
    if length(DI_index) > 1
        DI_index = DI_index(1);
    end
    DI(1,n) = (DI_index - (shift/bin)/2 - 1) * bin / unique_condition_num(3);

    % Minus with Plus
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_MP;
    for q = 1:2*(shift/bin)
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
        yi_temp = yi_shift;
        for r = 2:length(yi_shift{3,n})
            yi_shift{3,n}(r) = yi_temp{3,n}(r-1);
        end
        yi_shift{3,n}(1) = yi_temp{3,n}(end);
    end
    cov_MP = cov_MP((shift/bin)/2 + 1 : length(xi) - (shift/bin)/2); % limit possible DI's to only plausible values (+/- 90 deg)
    DI_index = find(cov_MP == max(cov_MP));
    if length(DI_index) > 1
        DI_index = DI_index(1);
    end
    DI(2,n) = (DI_index - (shift/bin)/2 - 1) * bin / (unique_condition_num(3)*2);
           
    % Zero with Plus
    yi_shift = yi_shift_init;  % re-initialize
    clear cov_ZP;
    for q = 1:2*(shift/bin)
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
        yi_temp = yi_shift;
        for r = 2:length(yi_shift{3,n})
            yi_shift{3,n}(r) = yi_temp{3,n}(r-1);
        end
        yi_shift{3,n}(1) = yi_temp{3,n}(end);
    end
    cov_ZP = cov_ZP((shift/bin)/2 + 1 : length(xi) - (shift/bin)/2); % limit possible DI's to only plausible values (+/- 90 deg)
    DI_index = find(cov_ZP == max(cov_ZP));
    if length(DI_index) > 1
        DI_index = DI_index(1);
    end
    DI(3,n) = (DI_index - (shift/bin)/2 - 1) * bin / unique_condition_num(3);

end

return;