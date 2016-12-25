% SilverBimodalTest for circular data only, Aki, 09/28/2006, modified by GY
clear all
load -ascii ModeTestData.txt
i=find(ModeTestData<0);
ModeTestData(i)=ModeTestData(i)+360;
Nb=1000;
modelist=1:3;
[numModes, p]=modalityTestForYong(ModeTestData',modelist,Nb)
