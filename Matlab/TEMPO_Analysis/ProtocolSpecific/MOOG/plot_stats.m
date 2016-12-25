function [r_x, chi_x] = plot_stats(statfiledir, filename1)%, filename2)

f1 = fullfile(statfiledir, filename1);
fid = fopen(f1);
line = fgetl(fid);
cnt=1;
while (line~=-1)
    
     % format for batch files
    % PATH  FILE    
    
    
        
              
        spaces = isspace(line);
        space_index = find(spaces);
        
        %stats
        L=length(line);
        r_x(cnt) = sscanf(line(space_index(1)+ 1:space_index(2) - 1), '%f');
        chi_x(cnt) = sscanf(line(space_index(2):L),'%f');
        

        % go through the backdoor
        cnt = cnt + 1;
        line = fgetl(fid);
        line = fgetl(fid);
end
fclose(fid);

% f2 = fullfile(statfiledir, filename2);
% fid = fopen(f2);
% line = fgetl(fid);
% cnt=1;
% while (line~=-1)
%     
%      % format for batch files
%     % PATH  FILE    
%     
%     
%         
%               
%         spaces = isspace(line);
%         space_index = find(spaces);
%         
%         %stats
%         L=length(line);
%         r_y(cnt) = sscanf(line(space_index(1)+ 1:space_index(2) - 1), '%f');
%         chi_y(cnt) = sscanf(line(space_index(2):L),'%f');
%         
% 
%         % go through the backdoor
%         cnt = cnt + 1;
%         line = fgetl(fid);
%         line = fgetl(fid);
%     end
% fclose(fid);
select=isnan(chi_x);
chi_x(select)=0;
% chi_y(select)=0;
% select=isnan(chi_y);
% chi_x(select)=0;
% chi_y(select)=0;
% figure;
% hold on;
% plot(r_x,r_y,'+');
% plot(chi_x,chi_y,'rx');
% plot(0:.01:1,0:.01:1,'k');
% xlabel('Rectified Cosine Model');
% ylabel('Gaussian Model');
% legend('R^2','Chi^2', 0);