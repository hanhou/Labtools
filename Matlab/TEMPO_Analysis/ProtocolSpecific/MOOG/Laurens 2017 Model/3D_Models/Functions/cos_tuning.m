function r1 = cos_tuning(param, coord,u_azi)
    if nargin == 2
        u_ele = coord(1:5);
        u_azi = coord(6:13);
    else
        u_ele = coord ;
    end
    n1 = param(1); % nonlinearty
    azi1 = param(2); % preferred azi
    ele1 = param(3); % preferred ele

    %Squeezing Nonlinearity
    F = @(x,n)((exp(n.*x)-1)./n);
    
    %Surface grid
    [ele, azi] = meshgrid(u_ele, u_azi);  
    [x, y, z] = sph2cart(azi, ele, ones(size(azi)));
    
    %Cosine tuning
    [x1, y1, z1] = sph2cart(azi1, ele1, 1);
    if n1 == 0
        r1 = [x(:) y(:) z(:)]*[x1 y1 z1]' ;
    else
        r1 = F([x(:) y(:) z(:)]*[x1 y1 z1]', n1);
    end
    ma = F(1, n1);mi = F(-1, n1);
    r1 = r1 - (ma + mi)/2;
    r1 = 2*r1/(ma - mi);
%     disp([n1]);disp([mi ma]);disp(r1(1))
end