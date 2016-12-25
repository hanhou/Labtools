numfig=round(dim(1)/2);
clear m;

for m=1:1:10
figure(m+1);
% set(gca,'Position', [5,5 1000,680], 'Name', 'Envelope');    

i=m*2-1;

  subplot(4,2,1);
    plot(res_x_up{i}(:,:)','r.');
    hold on;
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title1);
    
    plot(res_y_up{i}(:,:)','b.');
    hold off;

    text (10,0.7,['Cell No.', num2str(i)]);
    
  subplot(4,2,3);
    plot(res_x_down{i}(:,:)','r.');
    hold on;    
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title2);
      
    plot(res_y_down{i}(:,:)','b.');
    hold off;    
 
    
    
  subplot(4,2,5);
    plot(res_x_left{i}(:,:)','r.');
    hold on;   
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title3);
    
    plot(res_y_left{i}(:,:)','b.');
    hold off;   
  
    
    
  subplot(4,2,7);
    plot(res_x_right{i}(:,:)','r.');
    hold on;  
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title4);
    
    plot(res_y_right{i}(:,:)','b.');
    hold off;  

i=[];
i=m*2;

  subplot(4,2,2);
    plot(res_x_up{i}(:,:)','r.');
    hold on;
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title1);
    
    plot(res_y_up{i}(:,:)','b.');
    hold off;

     text (10,0.7,['Cell No.', num2str(i)]);
    
  subplot(4,2,4);
    plot(res_x_down{i}(:,:)','r.');
    hold on;    
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title2);
      
    plot(res_y_down{i}(:,:)','b.');
    hold off;    
 
    
    
  subplot(4,2,6);
    plot(res_x_left{i}(:,:)','r.');
    hold on;   
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title3);
    
    plot(res_y_left{i}(:,:)','b.');
    hold off;   
  
    
    
  subplot(4,2,8);
    plot(res_x_right{i}(:,:)','r.');
    hold on;  
    xlim( [1, 400] );
    set(gca, 'XTickMode','manual');
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
    ylim([-1,1]);
    ylabel('(deg)');
    title(title4);
    
    plot(res_y_right{i}(:,:)','b.');
    hold off;   

end