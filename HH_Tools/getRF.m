%% HH20130904

function getRF

if ishandle(77); close(77);end;   figure(77);
set(gca,'xtickmode','a');
axis equal; box off; axis([-45 45 -45 45]);

S.r1 = rectangle('position',[-45,-45,90,90]); hold on; 
S.l1 = line([-45 45],[0 0]);
S.l2 = line([0 0],[-45 45]); 
data.point = [];
data.rect_h = [];
set(S.r1,'userdata',data);
set([gca S.l1 S.l2 S.r1],'ButtonDownFcn',{@addPoint,S});  % Pass the handles to callbacks

function addPoint(~,~,S)

pos = get(gca,'CurrentPoint');
data = get(S.r1,'userdata'); % We attach userdata to r1

if (isempty(data.point)) 
    if ~isempty(data.rect_h)
        delete([data.rect_h data.point_h]);
        data.rect_h = [];
    end
    
    data.point = [pos(1,1) pos(1,2)];
    data.point_h(1) = plot(pos(1,1), pos(1,2),'.k');    
    set(gcf,'WindowButtonMotionFcn',{@moving,S});

elseif (data.w >0 && data.h >0)
    data.point_h(2) = plot(pos(1,1), pos(1,2),'.k');
    s = sprintf('[%2.0f, %2.0f, %2.0f, %2.0f]',data.x+data.w/2,data.y+data.h/2,data.w,data.h);
    clipboard('copy',s);
    fprintf('%s\n',s);
    data.point = [];
    set(gcf,'WindowButtonMotionFcn',{});
end

set(S.r1,'userdata',data);   % Pass back!!


function moving(~,~,S)

data = get(S.r1,'userdata');

if ~isempty(data.point)

    if ~isempty(data.rect_h)
        delete(data.rect_h); 
        data.rect_h = [];
    end;
    
    pos = get(gca,'CurrentPoint');
    data.x = min(pos(1,1),data.point(1));
    data.y = min(pos(1,2),data.point(2));
    data.w = abs(pos(1,1) - data.point(1));
    data.h = abs(pos(1,2)-data.point(2));
    
    if (data.w > 0 && data.h > 0)
        data.rect_h = rectangle('position',[data.x data.y data.w data.h]);
        set(data.rect_h,'ButtonDownFcn',{@addPoint,S});  % This is important!!
    end
    
    set(S.r1,'userdata',data); % Pass back!!
end
   