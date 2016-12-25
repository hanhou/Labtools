function TransRelativeDisparity_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	switch(Analysis{1})
	   case 'Plot Tuning Curves'
         TransRelativeDisparityCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
	   case 'Plot Vergence Angle'
         PlotTransVergenceAngle(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Plot Rasters/Histograms'
         PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, Protocol);
         disp('plot rasters here');
     case 'Fit Gabor'
         GaborFit_TransRelDisp(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
   end

return;