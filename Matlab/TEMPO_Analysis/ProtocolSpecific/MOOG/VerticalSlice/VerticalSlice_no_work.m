function VerticalSlice(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE) 
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

%  vertical slice code, begins here.  this code will calculate the vector sum in a 
%  in a single plane and then will find the difference in vector direction
%  and finally plot the distribution of vector directions for a given population
%  of cells, enjoy. JRB 4/16/07

vest = 1;
vis = 2;
combined =3;

%  possibility is for the different points that you want to look at, change
%  this if you are looking at a different plane, this is for the vertical
%  (fronto-parallel plane)
possibility = [0 45 90 135 180 225 270 315; 0 0 0 0 0 0 0 0 ];


vestibular_index = find(stim_type == vest & amplitude ~= 0 & (azimuth == 0 | azimuth == 180 | azimuth == 45 ...
                      | azimuth == 90 | azimuth == 135 | azimuth == 225 | azimuth == 270 | azimuth == 315 ...
                      ) & (elevation == 0));

visual_index = find(stim_type == vis & amplitude == 0 & (azimuth == 0 | azimuth == 180 | azimuth == 45 ...
                      | azimuth == 90 | azimuth == 135 | azimuth == 225 | azimuth == 270 | azimuth == 315 ...
                      ) & (elevation == 0));

% i like to have all of my data in one place thus i combine all the data
combined_matrix = cat(1,azimuth,elevation,spike_rates,stim_type,amplitude);

% this decomposes the above matrix into each of the conditions that are
% being compared using the above indices
% for i=1:length(combined_index)
%     combined_results(1:3,i) = combined_matrix(1:3,combined_index(i));
% end

for i=1:length(visual_index)
    visual_results(1:3,i) = combined_matrix(1:3,visual_index(i));
end

for i=1:length(vestibular_index)
    vestibular_results(1:3,i) = combined_matrix(1:3,vestibular_index(i));
end

% among all the visual trials this finds the ones that match the correct
% directions and writes them in index_sim such that every row responds to a
% column in possibility and every time it finds a match it writes it in the
% appropriate column
for j=1:length(visual_results)
    for i=1:length(possibility)
        if (visual_results(1:2,j) == possibility(:,i))
            index_sim(i,j) = j;
        end
    end
end

% makes a matrix for later anova testing using the smallest direction
% repition
for i=1:length(possibility)
    dimensions(i) = length(find(index_sim(i,:)));
    final_dim = min(dimensions);
end
anova_matrix_vis = zeros(length(possibility),final_dim);
[m,n] = size(anova_matrix_vis);


%calculates the mean firing rate for repeated trials
for i=1:length(possibility)  
    inter_vis = find(index_sim(i,:)~=0);
    average_vis(1:2,i) = visual_results(1:2,inter_vis(1));
    for j=1:length(inter_vis)
        fr_s(j) = visual_results(3,inter_vis(j));
        anova_matrix_vis(i,j) = visual_results(3,inter_vis(j));
    end
    average_vis(3,i) = mean(fr_s);
end
anova_matrix_vis = anova_matrix_vis';
p_vis = anova1(anova_matrix_vis,[],'off');

% this is slightly counterintuitive.  in order to specify each direction in
% the vertical plane, i summed each point, and since they are unique i used
% the sum to specify their direction in one plane, for example:  (azi,ele)
% for a motion that simulates straight up (0,-90) i summed x and y to equal
% -90, and in the plane that i was plotting this represented 90 degrees.
% this will obviously have to change to your convention
for i=1:length(average_vis)
    directions(i) = average_vis(1,i) + average_vis(2,i);
end

for i=1:length(directions)
    switch directions(i)
        case 0
            direction(i) = 0;
        case 45
            direction(i) = 45;
        case 90
            direction(i) = 90;
        case 135
            direction(i) = 135;
        case 180
            direction(i) = 180;
        case 225
            direction(i) = 225;
        case 270
            direction(i) = 270;
        case 315
            direction(i) = 315;
    end
end

% makes and input matrix with averaged firing rates and the resulting
% directions
input4vectorsum_vis = cat(1,direction,average_vis(3,:));

% converts from degrees to radians
for i=1:length(input4vectorsum_vis)
    input4vectorsum_conv_vis(i) = input4vectorsum_vis(1,i) * pi/180;
end

input4vectorsum_rad_vis = cat(1,input4vectorsum_conv_vis,average_vis(3,:));

% the actual vector summation begins here
for i=1:length(input4vectorsum_rad_vis)
    [x(i),y(i)] = pol2cart(input4vectorsum_rad_vis(1,i),input4vectorsum_rad_vis(2,i));
end

xy_vis = cat(1,x,y);

sum1 = sum(xy_vis,2);

[p,q] = cart2pol(sum1(1,1),sum1(2,1));

vectorsum_vis = cat(1,p,q);

% the following code is the same as above just specified for vestibular
% conditions

for j=1:length(vestibular_results)
    for i=1:length(possibility)
        if (vestibular_results(1:2,j) == possibility(:,i))
            index_sim_vest(i,j) = j;
        end
    end
end

for i=1:length(possibility)
    dimensions(i) = length(find(index_sim_vest(i,:)));
    final_dim = min(dimensions);
end
anova_matrix_vest = zeros(length(possibility),final_dim);
[m,n] = size(anova_matrix_vest);

for i=1:length(possibility)  
    inter_vest = find(index_sim_vest(i,:)~=0);
    average_vest(1:2,i) = vestibular_results(1:2,inter_vest(1));
    for j=1:length(inter_vest)
        fr_s_vest(j) = vestibular_results(3,inter_vest(j));
        anova_matrix_vest(i,j) = vestibular_results(3,inter_vest(j));
    end
    average_vest(3,i) = mean(fr_s_vest);
end

anova_matrix_vest = anova_matrix_vest';
p_vest = anova1(anova_matrix_vest,[],'off');

for i=1:length(average_vest)
    directions(i) = average_vest(1,i) + average_vest(2,i);
end

for i=1:length(directions)
    switch directions(i)
        case 0
            direction(i) = 0;
        case 45
            direction(i) = 45;
        case 90
            direction(i) = 90;
        case 135
            direction(i) = 135;
        case 180
            direction(i) = 180;
        case 225
            direction(i) = 225;
        case 270
            direction(i) = 270;
        case 315
            direction(i) = 315;
    end
end

input4vectorsum_vest = cat(1,direction,average_vest(3,:));

for i=1:length(input4vectorsum_vest)
    input4vectorsum_conv_vest(i) = input4vectorsum_vest(1,i) * pi/180;
end

input4vectorsum_rad_vest = cat(1,input4vectorsum_conv_vest,average_vest(3,:));

for i=1:length(input4vectorsum_rad_vest)
    [x(i),y(i)] = pol2cart(input4vectorsum_rad_vest(1,i),input4vectorsum_rad_vest(2,i));
end

xy_vest = cat(1,x,y);

sum1 = sum(xy_vest,2);

[p,q] = cart2pol(sum1(1,1),sum1(2,1));

vectorsum_vest = cat(1,p,q);

figure(2)
subplot(3,1,1),
polar([vectorsum_vis(1,1) vectorsum_vis(1,1)], [0 vectorsum_vis(2,1)],'r-')
hold on
polar([input4vectorsum_rad_vis(1,1) input4vectorsum_rad_vis(1,1)], [0 input4vectorsum_rad_vis(2,1)],'b-')
polar([input4vectorsum_rad_vis(1,2) input4vectorsum_rad_vis(1,2)], [0 input4vectorsum_rad_vis(2,2)],'b-')
polar([input4vectorsum_rad_vis(1,3) input4vectorsum_rad_vis(1,3)], [0 input4vectorsum_rad_vis(2,3)],'b-')
polar([input4vectorsum_rad_vis(1,4) input4vectorsum_rad_vis(1,4)], [0 input4vectorsum_rad_vis(2,4)],'b-')
polar([input4vectorsum_rad_vis(1,5) input4vectorsum_rad_vis(1,5)], [0 input4vectorsum_rad_vis(2,5)],'b-')
polar([input4vectorsum_rad_vis(1,6) input4vectorsum_rad_vis(1,6)], [0 input4vectorsum_rad_vis(2,6)],'b-')
polar([input4vectorsum_rad_vis(1,7) input4vectorsum_rad_vis(1,7)], [0 input4vectorsum_rad_vis(2,7)],'b-')
polar([input4vectorsum_rad_vis(1,8) input4vectorsum_rad_vis(1,8)], [0 input4vectorsum_rad_vis(2,8)],'b-')
title('Visual Response')
hold off

subplot(3,1,2),
polar([vectorsum_vest(1,1) vectorsum_vest(1,1)], [0 vectorsum_vest(2,1)],'r-')
hold on
polar([input4vectorsum_rad_vest(1,1) input4vectorsum_rad_vest(1,1)], [0 input4vectorsum_rad_vest(2,1)],'b-')
polar([input4vectorsum_rad_vest(1,2) input4vectorsum_rad_vest(1,2)], [0 input4vectorsum_rad_vest(2,2)],'b-')
polar([input4vectorsum_rad_vest(1,3) input4vectorsum_rad_vest(1,3)], [0 input4vectorsum_rad_vest(2,3)],'b-')
polar([input4vectorsum_rad_vest(1,4) input4vectorsum_rad_vest(1,4)], [0 input4vectorsum_rad_vest(2,4)],'b-')
polar([input4vectorsum_rad_vest(1,5) input4vectorsum_rad_vest(1,5)], [0 input4vectorsum_rad_vest(2,5)],'b-')
polar([input4vectorsum_rad_vest(1,6) input4vectorsum_rad_vest(1,6)], [0 input4vectorsum_rad_vest(2,6)],'b-')
polar([input4vectorsum_rad_vest(1,7) input4vectorsum_rad_vest(1,7)], [0 input4vectorsum_rad_vest(2,7)],'b-')
polar([input4vectorsum_rad_vest(1,8) input4vectorsum_rad_vest(1,8)], [0 input4vectorsum_rad_vest(2,8)],'b-')
title('Vestibular Response')
hold off

subplot(3,1,3),
polar([vectorsum_vis(1,1) vectorsum_vis(1,1)], [0 vectorsum_vis(2,1)],'r-')
hold on
polar([vectorsum_vest(1,1) vectorsum_vest(1,1)], [0 vectorsum_vest(2,1)],'b-')
title('Visual(red) and Vestibular(blue) Response')
hold off

saveas(gcf,['C:\MATLAB6p5\work\figures\' FILE '_slice.fig']);
figure(2);
close;

if ((p_vis > 0.05) | (p_vest > 0.05))
    difference = 9999;
else
    difference = 180/pi*(vectorsum_vis(1,1) - vectorsum_vest(1,1));
end

foldername = ('C:\MATLAB6p5\work\');
outfile1 = [foldername 'vertical_slice.dat'];

sprint_txt = ['%s\t %f\t %f\t %f\t'];
buff = sprintf(sprint_txt, FILE, difference, vectorsum_vis(1,1)*180/pi, vectorsum_vest(1,1)*180/pi);  


printflag = 0;
if (exist(outfile1, 'file') == 0)    %file does not yet exist
    printflag = 1;
end

fid = fopen(outfile1, 'a');

fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
return