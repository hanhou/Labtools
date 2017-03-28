 function norm_pos = GetPos()
% Use RBBOX to establish a position for a GUI component.
% object_handle is a handle to a uicomponent that uses
% any Units. Internally, figure Units are used.

h = gcf;
init_units = get(h,'Units'); % Cache initial position

set(h,'Units','pix');
f_pix_pos = get(gcf,'Position');

set(h,'Units','norm');
f_norm_pos = get(h,'Position');

disp(['=== Drag out a Position ===']);

waitforbuttonpress  % So that rbbox does not return immediately
norm_pos = rbbox;     % User drags out a rectangle, releases button
% Pressing a key aborts rbbox, so check for null width & height
norm_pos = fix(norm_pos * 1000)/1000;

if norm_pos(3) ~= 0 && norm_pos(4) ~= 0
    % Save and restore original units for object
    pix_pos = norm_pos .* [f_pix_pos(3) f_pix_pos(4) f_pix_pos(3) f_pix_pos(4)];
    fprintf('-- Object --\nNorm (clipboard): %s\nPix: %s\n',num2str(norm_pos),num2str(pix_pos));
    clipboard('copy', norm_pos);         % Place set string on system
else
%     set(object_handle,'Units', myunits);
    fprintf('-- Figure --\nNorm (clipboard): %s\nPix: %s\n',num2str(f_norm_pos),num2str(f_pix_pos));
    clipboard('copy', sprintf('set(gcf,''uni'',''norm'',''pos'',[%s]);',num2str(fix(f_norm_pos*1000)/1000)));         % Place set string on system
end

set(h,'Units', init_units);   % Restore initial units
