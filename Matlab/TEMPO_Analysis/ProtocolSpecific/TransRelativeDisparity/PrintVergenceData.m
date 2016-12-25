%---------------------------------------------------------------------------------------------------------------------
%-- PrintVergenceData.m : This function prints out Vergence Angle data for each surround disparity on the figure plot.
%--	KKW, 6/27/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintTransVergenceData(p_value, max_stats, min_stats, disp_groups, surr_disps, PATH, FILE);
	axis([0 100 0 100]);
	axis('off');
   font_size = 8;
   bump_size = 7;
   
   
   % type out stats onto screen
		xpos = -10;
      ypos = 20; 
      
      line = sprintf('SurrDisp        Pval');
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;

	for row = 1:size(surr_disps,1)
  		line = sprintf('%7.4f       %7.4f', surr_disps(row), p_value(row));
	   text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
   end
   
return;