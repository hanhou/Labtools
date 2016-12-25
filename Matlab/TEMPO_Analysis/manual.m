==========================
Manual of Combine Database
==========================

Introduction
============
This Matlab program will help you combine .log files and .htb(database) files. 
After combination, the log file and htb file still have one header.

Contents
========
1) Using selected Input Data Files to combine Database
2) Using Batch File to combine more than one Database
3) Rule of read Batch File
4) Rule of Output File Name generation
5) Select 'Skip FIXED_SEED checking'
6) Select 'Select Epoch(or Trial)'
7) Content in CombineDatabase.log file
8) Stop combining database
9) Select 'Show good trial number'
10) Manage database


1) Using selected Input Data Files to combine Database:
=======================================================
a)	Click 'Add Data File' button to add files that support Multi-Selection.
	If you add more than two files, output file will setup automatically.
b)	Click ‘Remove File’ button to remove files that also support Multi-Selection.
c)	Select ‘Output File Path’ by 'Browse' button
d)	Click 'Combine' button
e)	i) If we find anything different in header of input log files,
	then a window will pop up and show what different.
	ii) By default, we will use first input log file header contents,
	but you can modify in pop up window, so that we will use it in new combined log file.
	iii) If user click 'Cancel' in pop up window, we will not combine database.

For example,
Add two file:
	C:\Temp\m3c1203r13.htb 
	C:\Temp\m3c1203r14.htb 
Output File Path: C:\Temp
Output File Name: m3c1203r13_14

you get three files: m3c1203r13_14.log, m3c1203r13_14.htb and CombineDatabase.log


2) Using Batch File to combine more than one Database:
======================================================
a)	Select 'Input Data Files from Batch File' by 'Browse' button
b)	Select ‘Output File Path’ by 'Browse' button
c)	Click 'Combine' button

Example of Batch File:
...
%%% Post same day
Z:\Data\Tempo\Baskin\Raw\	m4c36r6.htb	'Behavioral Bootstrap'	...
%%% Post 1 days
Z:\Data\Tempo\Baskin\Raw\	m4c36r7.htb	'Behavioral Bootstrap'	...
Z:\Data\Tempo\Baskin\Raw\	m4c36r8.htb	'Behavioral Bootstrap'	...
%%% Post 2 days
Z:\Data\Tempo\Baskin\Raw\	m4c36r11.htb 	'Behavioral Bootstrap'	...
...


3) Rule of read Batch File:
===========================
a)	We will check each line first character.
b)	If there is more than or equal to two continue lineS that first character is 'Z' or 'z', 
	than we will combine all those files.
c)	Output File Name will assign automatically.


4) Rule of generation of Output File Name:
==========================================
a)	Each file name should include character 'm', 'c' and 'r' with number after them. i.e. m3c1203r13
b)	output file name still includes character 'm', 'c' and 'r' with number after them.
c)	If number is different from input files after character 'm', 'c' and 'r', 
	then combine number with '_'.

For example,
Input files: m3c1203r13
	     m4c1204r14 
Output file: m3_4c1203_1204r13_14 	


5) Select 'Skip FIXED_SEED checking':
=====================================
Since each log file has different FIXED_SEED, you can skip to check that and combine file quickly.


6) Select 'Select Epoch(or Trial)':
===================================
If you select the option of 'Select Epoch(or Trial)', then a pop up window will ask you input 
Start Epoch number and End Epoch number, that equal to trial# in log file.
It will only combine your selected Epoch number according to different input files.
It also will take more time to finish combining database.


7) Content of CombineDatabase.log file:
=======================================
a)	This file will exist in 'Output File Path'.
b)	Recording date of creating combined database file.
c)	Recording name and path of input files.
d)	Recording name and path of Output files.
e)	Recording different between input log files headers.
f)	Recording User selected new header parameters.
g)	Recording selected Epoch numbers, if user checked 'Select Epoch(or Trial)'
h)	If user cancels to combine database, it still has a record.
i)	We will append new record to the end of the file and will not discard existing contents.
Notice: You should open CombineDatabase.log file by WordPad or Word application.


8) Stop combining database:
===========================
When you click 'Stop' button, it will not stop immediately,
because Matlab will accumulate the commands before the stop command.  

9) Select 'Show good trial number':
===================================
When you selected 'Select Epoch(or Trial)', you also can select 'Show good trial number'.
It will show a good trial information in a new window and help you to select trial number.

10) Manage database:
===================
When you click 'Manage Data' button, it will create a new log and htb file for each input database,
that according to your selected trial number. That means you can resize the database.
Notice: Make sure you select different 'Output File Path' and Don't overwrite your original files.
        If not, you will lose the original files.
