datafile = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\HDispTuning\MUHDispCurves.dat';
fid = fopen(datafile, 'r');

pedestals = [0.1 0.25 0.025 -0.225 0.1 -0.1 0.05 0.15 -0.05 0.15 0.15 0.1 -0.05 0.2 0.15 0.15 0.15 0.2 0.3 0 0.2 0.15 0.05 0 -0.2 0.1 0.3 0 -0.1 0.2 0.3 0.1 0 0.05 0.1 0.3 0.2 -0.4 -0.1 0.2 0.2 0.2 0.15 0.2 0.15 0.15 0.2];

header = fgetl(fid);
index = 1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %disp(tline)
    [name, num, datastr] = strread(tline, '%s %f %[^\n]', 'bufsize', 10000);
    fname{index} = name;
    intnum = round(num);
    Npts(index) = intnum;

    all_data = sscanf(datastr{1},'%f');
    
    hdisp(index,:) = all_data(1:intnum)';
    mean_resp(index,:) = all_data(intnum+1:2*intnum)';
    err(index,:) = all_data(2*intnum+1:3*intnum)';
   
    index = index + 1;
end
fclose(fid);

%normalize the mean response data and errors
mm = max(mean_resp')';  %max response for each curve
divterm = (mm*ones(1,size(mean_resp,2)));
norm_mean_resp = mean_resp./divterm;

% figure;
% splot_ind=1;
% for i=1:3:(index-1)
%     subplot(8,6,splot_ind);
%     plot(hdisp(i:i+2,:)', mean_resp(i:i+2,:)', '-')
%     splot_ind = splot_ind + 1;
% end

figure;
splot_ind=1;
for i=1:3:(index-1)
    subplot(8,6,splot_ind);
    plot(hdisp(i:i+2,:)', norm_mean_resp(i:i+2,:)', '-');
    min_val =  min(min(norm_mean_resp(i:i+2,:)));
    hold on;
    plot([pedestals(splot_ind) pedestals(splot_ind)], [0 1], 'k-');
    ylim([0.9*min_val 1]);
    %ylim([0.2 1]);
    hold off;
    %axis('off');
    set(gca,'XTickLabel',{})
    set(gca,'YTickLabel',{})
    splot_ind = splot_ind + 1;
end

%plot(hdisp(1:3,:)', mean_resp(1:3,:)', '-')

%mm = max(mean_resp')'
%divterm = (mm*ones(1,9))
%norm_mean_resp = mean_resp./divterm
