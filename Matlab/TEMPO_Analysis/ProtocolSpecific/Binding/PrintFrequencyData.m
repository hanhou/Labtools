%---------------------------------------------------------------------------------------------------------------------
%-- PrintFrequencyData.m : THis function prints out spatial frequency or temporal frequency tuning curve data on the 
%--     current figure plot.
%--	BJP, 7/26/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintFrequencyData(p_value, base_rate, spont_level, amplitude, pref_freq, max_rate, width)
	axis([0 100 0 100]);
	axis('off');
   font_size = 9;
   bump_size = 7;
	
   xpos = -10;   
   ypos = 20;
   
   line = sprintf('Fitted Tuning Curve:');
	text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
   line = sprintf('Preferred Frequency = %1.5g   Resp = %0.5g', pref_freq, max_rate);
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
   line = sprintf('Base Rate = %0.5g   Spont = %0.5g', base_rate, spont_level);
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('Amplitude = %0.3g', amplitude);
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
   line = sprintf('Width (FWHM) = %0.3g', width);
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
   line = sprintf('ANOVA: P = %0.3g', p_value);
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
  
 return;