%---------------------------------------------------------------------------------------------------------------------
%-- PrintGeneralData.m : THis function prints out some useful general information at the top of each figure plot.
%--	GCD, 1/26/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintCircTransData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Tempo_Defs;		%needed for defns
ProtocolDefs;	%references to protocol Names 1/4/01 BJP

axis([0 100 0 100]);
axis('off');
xpos = -10;
ypos = 200;
font_size = 9;
bump_size = 7;

	line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
	text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


return;