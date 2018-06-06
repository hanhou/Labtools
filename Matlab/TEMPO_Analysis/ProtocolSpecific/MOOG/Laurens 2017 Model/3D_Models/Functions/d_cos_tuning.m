function r = d_cos_tuning(param, coord)
    u_ele = coord(1:5);
    u_azi = coord(6:13);
    
    n1 = param(1);
    azi1 = param(2);
    ele1 = param(3);
    
    n2 = param(4);
    azi2 = param(5);
    ele2 = param(6);
    
    a = param(7);
    c = param(8);
    
    %Squeezing Nonlinearity
    F = @(x,n)((exp(n.*x)-1)./n);
    
    %Surface grid
    [ele, azi] = meshgrid(u_ele, u_azi);
    [x, y, z] = sph2cart(azi, ele, ones(size(azi)));
    
    %Cosine tuning 1
    [x1, y1, z1] = sph2cart(azi1, ele1, 1);
    r1 = F([x(:) y(:) z(:)]*[x1 y1 z1]', n1);
    r1 = r1 - (max(r1) + min(r1))/2;
    r1 = 2*r1/(max(r1) - min(r1));
    
    %Cosine tuning 2
    q_s = cos(azi2/2);
    q_v = sin(azi2/2)*[x1 y1 z1];
    q = [q_v q_s];
    
    [x2, y2, z2] = sph2cart(azi1, ele1+ele2, 1);
    v = qvqc(q, [x2 y2 z2]');
    x2 = v(1);
    y2 = v(2);
    z2 = v(3);
    
    r2 = F([x(:) y(:) z(:)]*[x2 y2 z2]', n2);
    r2 = r2 - (max(r2) + min(r2))/2;
    r2 = 2*r2/(max(r2) - min(r2));
    
    r = r1 + a*r2;
    r = r - (max(r) + min(r))/2;
    r = 2*r/(max(r) - min(r));
    r = r + c;
end