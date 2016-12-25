function SurroundMapping_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

	switch(Analysis{1})
	   case 'Plot Surround Response Map'
         SurroundMappingCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
         %case 'Gradient Predictor'
         %SurroundModel(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
      case 'Plot Rasters/Histograms'
   		disp('plot rasters here');   
      case 'Plot Vergence Data' 
   		disp('plot vergence here');   
   end

return;