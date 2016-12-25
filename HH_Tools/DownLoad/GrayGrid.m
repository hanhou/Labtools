function hAx2 = GrayGrid(hAx1)

hFig = get(hAx1,'parent');

%# create a second transparent axis, as a copy of the first
hAx2 = copyobj(hAx1,hFig);
delete( get(hAx2,'Children') )
set(hAx2, 'Color','none', 'Box','on', ...
    'XGrid','off', 'YGrid','off')

%# show grid-lines of first axis, style them as desired,
%# but hide its tick marks and axis labels
set(hAx1, 'XColor',[0.9 0.9 0.9], 'YColor',[0.9 0.9 0.9], ...
    'XMinorGrid','on', 'YMinorGrid','on', 'MinorGridLineStyle','-', ...
    'XTickLabel',[], 'YTickLabel',[]);
xlabel(hAx1, ''), ylabel(hAx1, ''), title(hAx1, '')

%# link the two axes to share the same limits on pan/zoom
linkaxes([hAx1 hAx2], 'xy');
