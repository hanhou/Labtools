%---------------------------------------------------------------------------------------------------------------------
%-- PrintRelDispData.m : THis function prints out H Disp tuning curve specific data on the igure plot.
%--	BJP, 2/1/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintRelDispData(data, p_value, avg_resp, max_stats, min_stats, disp_groups, spont_level, surr_disps, PATH, FILE);

	TEMPO_Defs;

	axis([0 100 0 100]);
	axis('off');
   font_size = 8;
   bump_size = 7;
   
   %get the value of the center patch direction in the dots_params matrix
	ctr_dir = data.dots_params(DOTS_DIREC,1,PATCH1);
   
	%get the value of the center patch speed in the dots_params matrix
   ctr_spd = data.dots_params(DOTS_SPEED,1,PATCH1);
   
   %get the value of the center patch coherence in the dots_params matrix
   ctr_coh = data.dots_params(DOTS_COHER,1,PATCH1);
   
   %get the value of the center patch size in the dots_params matrix
   ctr_siz = data.dots_params(DOTS_AP_XSIZ,1,PATCH1);
   
	%get the value of the surround patch direction in the dots_params matrix
   surr_dir = data.dots_params(DOTS_DIREC,1,PATCH4);
   
   %get the value of the surround patch speed in the dots_params matrix
   surr_spd = data.dots_params(DOTS_SPEED,1,PATCH4);
   
   %get the value of the surround patch coherence in the dots_params matrix
   surr_coh = data.dots_params(DOTS_COHER,1,PATCH4);
   
   %get the value of the surround patch size in the dots_params matrix
   surr_siz = data.dots_params(DOTS_AP_XSIZ,1,PATCH4);



   
   % type out stats onto screen
		xpos = -10;
      ypos = 20; 
      
      line = sprintf('SurrDisp        DTI         Pval');
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;

	for row = 1:size(surr_disps,1)
	   % first calculate some metrics and stats
   	DTI = 1 - (min_stats(row).y - spont_level)/(max_stats(row).y - spont_level);
           
   	
      
  		line = sprintf('%7.4f       %7.4f       %7.4f', surr_disps(row), DTI, p_value(row));
	   text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;

      
      end


      xpos = 60;
      ypos = 20;
      
      line = sprintf('CtrDir     CtrSpd      CtrCoh');
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
      
      line = sprintf('%7.4f   %7.4f   %7.4f', ctr_dir, ctr_spd, ctr_coh);
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
      
      line = sprintf('SurrDir    SurrSpd     SurrCoh');
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
      
      line = sprintf('%7.4f   %7.4f   %7.4f', surr_dir, surr_spd, surr_coh);
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
      
      line = sprintf('CtrSiz     SurrSiz');
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
      
      line = sprintf('%7.4f   %7.4f', ctr_siz, surr_siz);
      text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;




   
return;