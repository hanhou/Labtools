%----------------------------------------------------------------------------------------------------------
%-- TEMPO_Defs.m: This file contains a bunch of definitions that are needed to interpret and work with
%--	data from TEMPO.  Some of the indices here MUST match up with their respective indices in the TEMPO
%--	protocol.  I have tried to put all such TEMPO-specific defines in this module, and urge you to 
%--	maintain this structure.   GCD, 1/3/2000
%----------------------------------------------------------------------------------------------------------

%defines for indexing the different databases; these are potentially protocol-specific
%NOTE: these indices must match the order of the databases defined in the TEMPO protocol!
EYE_DB = 1;		%eye movement samples
SPIKE_DB = 2;	%spike times
EVENT_DB = 3;	%event times
LFP_DB = 4;		%LFP samples

%DEFINES RELATED TO DATABASES

%eye channels (Default settings, could be overriden. HH20150722)
LEYE_H = 1;	
LEYE_V = 2;
REYE_H = 3;
REYE_V = 4;

%D/A inputs for moving fixation window command signal
DA_H = 5;
DA_V = 6;

%spike channels
SPIKE1 = 1;
SPIKE2 = 2;
SPIKE3 = 3;
SPIKE4 = 4;
%event codes (these MUST MATCH the list in the TEMPO protocol, ecodes.PRO)
TRIAL_START_CD		= 1;		%trial started   01
FP_ON_CD		    = 2;	    %Fixation Point on   02
IN_FIX_WIN_CD		= 3;		%Entered fixation window  03
VSTIM_ON_CD 		= 4;		%visual stimulus on  04
VSTIM_OFF_CD 		= 5;		%visual stimulus off  05
TARGS_ON_CD			= 6;		%targets on  06
SACCADE_BEGIN_CD	= 7;		%saccade has begun; monkey left fixation window  07
IN_T1_WIN_CD		= 8;		%monkey is in the T1 window  08
IN_T2_WIN_CD		= 9;		%monkey is in the T2 window   09
BROKE_FIX_CD		= 10;		%monkey broke fixation     0A
BROKE_VERG_CD		= 11;		%monkey broke vergence     0B
SUCCESS_CD			= 12;		%trial was successful       0C
REWARD_CD			= 13;		%reward is delivered    0D
PUNISH_CD			= 14;		%beep is delivered    0E
TRIAL_END_CD		= 15;		%trial ended     0F
MICROSTIM_ON_CD 	= 16;		%microstim tuned on   10
MICROSTIM_OFF_CD 	= 17;		%microstim tuned off   11
CUE_ON_CD			= 18;       %cue is turned on   12
CUE_OFF_CD			= 19;       %cue is turned off   13
AXIS_CUE_ON_CD		= 20;       %axis cue is turned on   14
AXIS_CUE_OFF_CD		= 21;       %axis cue is turned off  15
SUB_TRIAL_END_CD    = 22;       %end of subtrial (in multistim protocols)   16
FP_OFF_CD           = 23;       %fixation point off and nothing else happens   17


%event description strings
%NOTE: the order and number of these strings has to match the list of event codes above
event_names = ...
   {  'Trial Start',
      'Fixation Point ON',
      'Entered Fixation Window',
      'Visual Stimulus ON',
      'Visual Stimulus OFF',
      'Targets ON',
      'Left Fixation Window',
      'Entered T1 Window',
      'Entered T2 Window',
      'Broke Fixation',
      'Broke Vergence',
      'Trial Successful',
      'Reward ON',
      'Punish (beep) ON',
      'Trial End',
      'Microstim ON',
      'Microstim OFF',
      'Cue ON',
      'Cue OFF',
      'Axis Cue ON',
      'Axis Cue OFF' 
      };

%these are possible values of the OUTCOME variable that is read in from the TEMPO log file
%(see misc_keywords below).  Note that this list of values *MUST* correspond to those in the 
%TEMPO paradigm.  GCD 1/6/00
CORRECT = 0;			%correct trial
ERR_ACQUIRE_FP = 1;		%set when fails to acquire fixation point
ERR_BROKE_FIXATION = 2;	%set when monkey breaks fixation during trial
ERR_NO_CHOICE = 3;		%value if monkey fails to leave FIxWin after targets shown
ERR_MISSED_TARGETS = 4;	%left FixWin but missed both target windows
ERR_WRONG_CHOICE = 5;	%made a saccade to a target, but wrong one
ERR_BROKE_VERGENCE = 6;	%set when monkey breaks vergence criterion during trial

ProtocolDefs;	%Call defns of keywords that were previously in this file BJP 1/4/01

% Different rigs have different rig's input variable define.
% Redefinition here will promise to give a correct setting of variable
RigHandle = findobj(gcbf, 'Tag', 'RigPopup');
if ~isempty(RigHandle) 
    rigger = get(RigHandle, 'Value');
    switch(rigger)
        case 1 % Moog 1
            DA_H = 7;
            DA_V = 8;
            ACC_X = 5;
            ACC_Y = 6;
        case 2 % Moog 2
            DA_H = 5;
            DA_V = 6;
            ACC_X = 7;
            ACC_Y = 8;
        case 3 % Moog 3
            DA_H = 7;
            DA_V = 8;
            ACC_X = 5;
            ACC_Y = 6;
        case 4 % Mini Moog
            DA_H = 7;
            DA_V = 8;
            ACC_X = 5;
            ACC_Y = 6;
    end
end


