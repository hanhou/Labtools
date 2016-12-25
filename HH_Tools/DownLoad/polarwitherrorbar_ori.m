function [] = polarwitherrorbar(angle,avg,error)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first two input variables ('angle' and 'avg') are same as the input 
% variables for a standard polar plot. The last input variable is the error
% value. Note that the length of the error-bar is twice the error value we
% feed to this function. 
% In order to make sure that the scale of the plot is big enough to
% accommodate all the error bars, i used a 'fake' polar plot and made it
% invisible. It is just a cheap trick. 
% The 'if loop' is for making sure that we dont have negative values  when
% an error value is substrated from its corresponding average value. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_data = length(angle);

fake = polar(angle,max(avg+error)*ones(size(angle))); set(fake,'Visible','off'); hold on; 

h_p = polar(angle,avg,'-b');
set(h_p,'LineWidth',3);

for ni = 1 : n_data
    if (avg(ni)-error(ni)) < 0
        h_p = polar(angle(ni)*ones(1,3),[0, avg(ni), avg(ni)+error(ni)],'-b'); 
    else
        h_p = polar(angle(ni)*ones(1,3),[avg(ni)-error(ni), avg(ni), avg(ni)+error(ni)],'-b'); 
    end
    set(h_p,'LineWidth',3);
end

hold off
