%% Giving a matrix and calculate its vector_sum , used for calculate
%% direction selectivity in 3D space  -GY

function [vector_azi, vector_ele, vector_amp] = vectorsum(resp);

% define initial data
pi=3.14159;
% initiates each vector's projection onto x,y,z axis respectively
xs=0;     
ys=0;     
zs=0;
% remove the single dimension
%resp = squeeze(resp2);
% calculate azimuth and elevation
dimension = size(resp);
azimuth = 0: 360/ dimension(2) : (360-360/dimension(2));

if dimension(1) == 1   % added for 1D data, CRF
    elevation = zeros(1,dimension(2));
else
    elevation = -90: 180/ (dimension(1)-1) : 90; % (original line)
end 

% number of vectors
if dimension(1) == 1   % added for 1D data, CRF
    vector_num = dimension(2);
else
    vector_num = dimension(2) * (dimension(1)-2) + 2;  % (original line)
end
    
for i=1:dimension(2)
    for j=1:dimension(1)        
        [xn,yn,zn] = sph2cart(azimuth(i)*pi/180,elevation(j)*pi/180, resp(j,i) );
        xs = xs + xn;
        ys = ys + yn;
        zs = zs + zn;    
    end
end

% go back to calculate the vector sum's azimuth, elevation and amplitude 
[temp_azi,temp_ele,vectorsum_amp] = cart2sph( xs,ys,zs );
vectorsum_azi = temp_azi * 180/ pi;
vectorsum_ele = temp_ele * 180/ pi;

% if vectorsum azimuth is negtive
if (vectorsum_azi < 0)
    vectorsum_azi = vectorsum_azi + 360;
end

% output result
vector_azi = vectorsum_azi;
vector_ele = vectorsum_ele;
vector_amp = vectorsum_amp;


return;