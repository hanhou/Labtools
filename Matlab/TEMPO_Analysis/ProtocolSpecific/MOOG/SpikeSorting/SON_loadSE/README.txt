********* Readme for LoadSE_SON ***************************************************

LoadSE_SON is an MClust Loading engine for CED SON files (*.smr).  It requires the SON library 
for Matlab written by Malcolm Lidierth, which is contained in this archive.  Check the contributed 
software page at http://www.ced.co.uk for updates to the SON library.
 
To use, place LoadSE_SON.m in the MClust "LoadingEngines" directory and put the SON library directory 
in your Matlab path.  

This version is for use only with single electrode recordings. It uses the first wavemark channel 
encountered, so if you have more than one wavemark channel in your SMR file make sure the channel 
you want imported is the lowest numbered.  This engine can be easily modified to use up to four 
channels. As is, it zero pads the three unused channels for compatability with MClust.  The engine
makes no assumptions about the number of samples per spike, so it should work with any number of 
points that MClust supports.  

Tested with SON library version 1.02, Spike2 version 4.19, and MClust 3.3.  

Shane Heiney <shane@vor.wustl.edu>
17 March 2004



********* Original Readme file from SON library ***********************************

These m-files use features introduced with version 5 of MATLAB and will not work 
with versions earlier than that.
They were developed with MATLAB version 6.1 release 12 so there could be glitches
with versions earlier than that.


Directories
To use the m-files unzip SON.ZIP and either
[1] copy them to a directory on your current MATLAB path
or better
[2] place them in a directory of their own and add this to your MATLAB path. To 
have the new directory included by default each time MATLAB starts edit the 
startup.m file which is probably in the matlab "..\toolbox\local" directory  (it 
might be called startupsav.m if it is not currently used).
A line such as:
path(path,'c:\matlab6p1\work\son');
adds the new directory to the default path. Use 'help path' in MATLAB for 
details


You will probably need to edit the SONTest.m file to point to the correct 
directory for the CED SONFix.exe file if you use the SONTest function. Edit the 
line that has:
SONFIX='c:\spike403\sonfix.exe';


Malcolm Lidierth
25/03/02

