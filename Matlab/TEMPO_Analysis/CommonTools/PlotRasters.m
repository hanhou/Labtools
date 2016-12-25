function PlotRasters(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE, Protocol);

	line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
   Tempo_Defs;
   ProtocolDefs;
       
    [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);
   num_modality = length(param);
 
figure
for trial = BegTrial:EndTrial
    spikes = [];
    spikes2 = [];
    spikes =  find(data.spike_data(SpikeChan, (StartEventBin(trial) + StartOffset):(StopEventBin(trial) + StopOffset), trial) == 1 );   
    plot (spikes/1000,trial*ones(1,length(spikes)),'k.');
    spikes2 = find(data.spike_data(SpikeChan2, (StartEventBin(trial) + StartOffset):(StopEventBin(trial) + StopOffset), trial) == 1 );   
    plot (spikes2/1000,(trial + 0.5)*ones(1,length(spikes2)),'r.');
    hold on;
end
xlabel('Time (ms)');
ylabel('Trial Number');


figure

yy = 0;
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:size(conditions,2);												% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
   %analyze by unique condition
	for cond = 1: num_conditions
     	% These next few lines check for trials belonging to a single condition
     	SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
     	for modality = 1: num_modality         	
            NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
      	    SetTrials = SetTrials & NextSetTrials;				
        end	 
     
        reps = find(SetTrials==1); 
        num_reps = length(reps);		   
      
        spikes = [];    
        incr = (1 - 0.2)/(num_reps*2)
 	     % Now do cross correlograms 
   	   for trial = 1:num_reps
            spikes = [];
            spikes2 = [];

            spikes1 = find(data.spike_data(SpikeChan, (StartEventBin(reps(trial) ) + StartOffset):(StopEventBin(reps(trial) ) + StopOffset), reps(trial) ) == 1 );   
            spikes2 = find(data.spike_data(SpikeChan2, (StartEventBin(reps(trial) ) + StartOffset):(StopEventBin(reps(trial) ) + StopOffset), reps(trial) ) == 1 );   
            plot (spikes1/1000,yy*ones(1,length(spikes1)),'k.');
            yy = yy + incr;
            plot (spikes2/1000,yy*ones(1,length(spikes2)),'r.');
            hold on;
            yy = yy + incr;
       end   
        yy = yy + 0.2;
   end
ylim([0,yy]);
xlabel('Time (ms)');
ylabel('Trials Sorted by Condition Number');

