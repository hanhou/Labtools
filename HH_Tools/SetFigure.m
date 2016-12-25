function SetFigure(size)
% Make figure more beautiful. @HH2013

if nargin == 0
    size = 15;
end

set(gcf,'color','w'); % Color
% set(findall(gcf,'-property','fontsize'),'fontsize',size); % Font size
set(findall(gcf,'fontsize',10),'fontsize',size); % Font size

set(findall(gcf,'tickdir','i'),'tickdir','o'); % Tick direction
set(findall(gcf,'type','axes'),'ticklength',[0.02 0],'LineWidth',1.5,'color','none'); % Tick length
	
% set(findall(gcf,'type','axes','linewidth',0.5,'-not','tag','legend'),'linewidth',15); % Tick width

set(findall(gcf,'type','axes','-not','tag','legend'),'box','off'); % Box off

