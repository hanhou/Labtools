function [] = polarwitherrorbar(angle,avg,error,marker,width,simple_mode)

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

% Add simple_mode. HH20180609

if nargin == 3
    marker = '-b';
    width = 2.5;
    simple_mode = 0;
elseif nargin == 4
    width = 2.5;
    simple_mode = 0;
elseif nargin == 5
    simple_mode = 0;
end

n_data = length(angle);

fake = polarInHeading(angle,max(avg+error)*ones(size(angle)),marker, simple_mode); set(fake,'Visible','off'); hold on; 

h_p = polarInHeading(angle,avg,marker,simple_mode );
set(h_p,'LineWidth',width);

for ni = 1 : n_data
    if (avg(ni)-error(ni)) < 0
        h_p = polarInHeading(angle(ni)*ones(1,3),[0, avg(ni), avg(ni)+error(ni)],marker,simple_mode ); 
    else
        h_p = polarInHeading(angle(ni)*ones(1,3),[avg(ni)-error(ni), avg(ni), avg(ni)+error(ni)],marker,simple_mode ); 
    end
    set(h_p,'LineWidth',width);
end

hold off
