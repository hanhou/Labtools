% Coherence analysis --YG, 07/30/2010
% %-----------------------------------------------------------------------------------------------------------------------
function ZCFUNC(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
TEMPO_Defs
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_a = data.moog_params(TRAPEZOID_ACC  ,:,MOOG);
temp_t1 = data.moog_params(TRAPEZOID_TIME  ,:,MOOG); 
temp_t2 = data.moog_params(TRAPEZOID_TIME3  ,:,MOOG); 
temp_v = temp_a.*temp_t1;
temp_d = temp_a.*temp_t1.*(temp_t1+temp_t2);
total_trials = length(temp_a);
LONG = 2;
SHORT = 1;
for i= 1 : total_trials
    temp = data.event_data(1,:,i + BegTrial-1);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = LONG;
    elseif (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = SHORT;
    else
     %   choice(i) = RIGHT;
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end
end
%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 

azimuth = temp_azimuth(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);

unique_a = munique(temp_a);
unique_t1 = munique(stimtemp_t1);
unique_t2 = munique(temp_t2);

spon_found = find(null_trials==1);     
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);


return;

