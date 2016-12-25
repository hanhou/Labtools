% DFT_boot.m:

function [p, p_vel, p_acc] = DFT_perm(time, data, dft_ratio, vel_dftr, acc_dftr, repeat, dftcutoff)

if sum(data) == 0 | dft_ratio == 0
    p = 1; p_vel = 1; p_acc = 1;
    return

else
	for i = 1:repeat
        shuf = randperm(length(data));
        shuf_dat = data(shuf);
        [f, amp, resp_phase] = FT(time, shuf_dat, length(time), 1, 0);
        f = round(f*100)/100; % get rid of some floating point issues
        f1 = mean(amp(find(f > 0 & f <= dftcutoff)));
        f2 = mean(amp(find(f > dftcutoff)));
        if f2 == 0
            dftr_shuf(i) = 0;
        else
            dftr_shuf(i) = f1/f2;
        end
        
        pol = -1;
        for s = find(f > 0 & f <= dftcutoff)
            amp_vel(s) = pol*amp(s)*cos(resp_phase(s));
            amp_acc(s) = pol*amp(s)*sin(resp_phase(s));
            pol = -pol;
        end
        v1 = mean(amp_vel(find(f > 0 & f <= dftcutoff)));
        if f2 == 0
            vel_dftr_shuf(i) = 0;
        else
            vel_dftr_shuf(i) = v1/f2;
        end
        a1 = mean(amp_acc(find(f > 0 & f <= dftcutoff)));
        if f2 == 0
            acc_dftr_shuf(i) = 0;
        else
            acc_dftr_shuf(i) = a1/f2;
        end
	end
	
	% p values
	p = 2*length(find(dftr_shuf >= dft_ratio))/repeat;
	if vel_dftr > 0
        p_vel = 2*length(find(vel_dftr_shuf >= vel_dftr))/repeat;
	else
        p_vel = 2*length(find(vel_dftr_shuf <= vel_dftr))/repeat;
	end
	if acc_dftr > 0
        p_acc = 2*length(find(acc_dftr_shuf >= acc_dftr))/repeat;
	else
        p_acc = 2*length(find(acc_dftr_shuf <= acc_dftr))/repeat;
	end

	return
    
end



% %%--------------MODIFIED DFT_PERM FOR BPR-ADJUSTED ---------------------%%
% %%--------------        April 27, 2008             ---------------------%%
% 
% %this takes in 2 sec of data and the original dft ratio of the data. It
% %permutes the data and recalculates DFT ratio.
% % MODIFIED BY CRF -- 9/2007: New vel- and acc-DFTR method, and some other fixes.
% function [p, p_vel, p_acc] = DFT_perm(time, data, dft_ratio, vel_dftr, acc_dftr, repeat, dftcutoff)
% 
% load('Z:\Data\Tempo\Batch Files\Suhrud\bpr_lookup_table.mat');
% 
% if sum(data) == 0 | dft_ratio == 0
%     p = 1; p_vel = 1; p_acc = 1;
%     return
%     
% else
% 	for i = 1:repeat
%         shuf = randperm(length(data));
%         shuf_dat = data(shuf);
%         [f, amp, resp_phase] = FT(time, shuf_dat, length(time), 1, 0);
%         f = round(f*100)/100; % get rid of some floating point issues
%         f1 = mean(amp(find(f > 0 & f <= dftcutoff)));
%         f2 = mean(amp(find(f > dftcutoff)));
%         if f2 == 0
%             dftr_shuf(i) = 0;
%         else
%             dftr_shuf(i) = f1/f2;
%         end
%         
%         pol = -1;
%         for s = find(f > 0 & f <= dftcutoff)
%             amp_vel(s) = pol*amp(s)*cos(resp_phase(s));
%             amp_acc(s) = pol*amp(s)*sin(resp_phase(s));
%             pol = -pol;
%         end
%         v1 = mean(amp_vel(find(f > 0 & f <= dftcutoff)));
%         if f2 == 0
%             vel_dftr_shuf(i) = 0;
%         else
%             vel_dftr_shuf(i) = v1/f2;
%         end
%         a1 = mean(amp_acc(find(f > 0 & f <= dftcutoff)));
%         if f2 == 0
%             acc_dftr_shuf(i) = 0;
%         else
%             acc_dftr_shuf(i) = a1/f2;
%         end
%         
%         dft_sum_shuf(i) = abs(vel_dftr_shuf(i)) + abs(acc_dftr_shuf(i));
%         
%         vel_pct_shuf(i) = vel_dftr_shuf(i) / dft_sum_shuf(i);
%         vel_pct_shuf(i) = round(vel_pct_shuf(i)*100) / 100;
%         acc_pct_shuf(i) = acc_dftr_shuf(i) / dft_sum_shuf(i);
%         acc_pct_shuf(i) = round(acc_pct_shuf(i)*100) / 100;
%         
%         % %  %%NEW STUFF APR 15th BEGIN  
%         datatogetmean = shuf_dat;
%         meanfiring = mean(shuf_dat(1:31));
%         sizeofbins = 10;
%         slider = 1;
%         startat = 1;
%         nn = 1;
%         endat = startat +sizeofbins-1;
%         while endat < length(shuf_dat)
%             if sum(datatogetmean(startat : endat) )~=0 
%                 meandata(nn) = mean(datatogetmean(startat : endat) );
%             else 
%                 meandata(nn) = 0;
%             end
%             startat = startat + slider;
%             endat = startat + sizeofbins -1;
%             nn = nn + 1;
%         end
%         maxmean = meandata(find(meandata == max(meandata)));
%         bpr_shuf(i) = meanfiring/maxmean(1);
%         bpr_shuf(i) = round(bpr_shuf(i)*100) / 100;
%         if bpr_shuf(i) >0.5
%             bpr_shuf(i) = 0.5;
%         end
% 
%         if isempty(find(bpr_lookup == bpr_shuf(i) & vel_pct_output == vel_pct_shuf(i)))
%             vel_pct_shuf_adj(i) = vel_pct_shuf(i);
%         else      
%             vel_pct_shuf_adj(i) = max(vel_pct_input( find(bpr_lookup == bpr_shuf(i) & vel_pct_output == vel_pct_shuf(i)) ));
%         end
%     
%         if isempty(find(bpr_lookup == bpr_shuf(i) & acc_pct_output == acc_pct_shuf(i)))
%             acc_pct_shuf_adj(i) = acc_pct_shuf(i);
%         else
%             acc_pct_shuf_adj(i) = max(acc_pct_input( find(bpr_lookup == bpr_shuf(i) & acc_pct_output == acc_pct_shuf(i)) ));
%         end
%         
%         vel_dftr_shuf(i) = vel_pct_shuf_adj(i) * dft_sum_shuf(i);
%         acc_dftr_shuf(i) = acc_pct_shuf_adj(i) * dft_sum_shuf(i);
% 
% 	end
% 
% 	% p values
% 	p = 2*length(find(dftr_shuf >= dft_ratio))/repeat;
% 
%     if vel_dftr > 0
%         p_vel = 2*length(find(vel_dftr_shuf >= vel_dftr))/repeat;
% 	else
%         p_vel = 2*length(find(vel_dftr_shuf <= vel_dftr))/repeat;
% 	end
% 	if acc_dftr > 0
%         p_acc = 2*length(find(acc_dftr_shuf >= acc_dftr))/repeat;
% 	else
%         p_acc = 2*length(find(acc_dftr_shuf <= acc_dftr))/repeat;
% 	end
% 
% 	return
%     
% end