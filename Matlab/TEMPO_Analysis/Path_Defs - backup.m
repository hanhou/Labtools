%----------------------------------------------------------------------------------------------------------
%-- Path_Defs.m: Path Definitions for various common and protocol-specific analysis routines.  You can
%--	simply modify BASE_PATH to change the base directory for all routines, e.g., when moving code
%--	to a different machine.  GCD, 1/3/2000
%----------------------------------------------------------------------------------------------------------

%Base path specification for protocol-specific analysis routines
BASE_PATH = 'Z:\LabTools\Matlab\TEMPO_Analysis\';

%add path for common tools
junk_str = [BASE_PATH 'CommonTools\'];
addpath(junk_str);

%add paths for individual protocol analysis tools
junk_str = [BASE_PATH 'ProtocolSpecific\DirectionTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SpeedTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\HDispTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\HGradTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\DepthDiscrim'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\DirectionDiscrim'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\AxisCuedDirectionTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\CuedDirectionDiscrim'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\OrientCueDirectionTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SizeTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SurroundMapping'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SurroundTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\RelativeDisparity'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\Stereoacuity'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\EyeCalibration'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\Binding'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\RFMapping'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\TransRelativeDisparity'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SimDistDispOnly'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SimDistVergOnly'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SimDistDispVerg'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SimDistCurvatureDiscrim'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\DepthDiscrim_NOVAR'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\DelayedSaccade'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MovingTarget'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning_Curvefit'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\DelayedSacc'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SurfOrientTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SurfSpdTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SurfDepthTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MotionParallaxFix'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\psignifit'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\TILT_TRANSLATION'];
addpath(junk_str);%added by AHC
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\Reverse Correlation'];
addpath(junk_str);%added by AHC
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\MU activity'];
addpath(junk_str);%added by AHC
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\VaryFixation'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination_2I'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\AzimuthTuning1D'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\Cueconflict2D'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\MemSacc']; %memory saccade
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\MOOG\Pursuit']; %visual pursuit
addpath(junk_str);