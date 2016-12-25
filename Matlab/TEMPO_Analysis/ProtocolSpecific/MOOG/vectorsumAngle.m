% new version, deals with uneven distributed azimuth or elevation 04/29/05 GY
% response is given by a vector resp, azimuth and elevation are also given
% by vector with azi and ele, respectively
% vector numbers should be the same

function [vector_azi, vector_ele, vector_amp] = vectorsumAngle(resp,azi,ele)

if nargin < 3   % Just azimuth
    ele = zeros(size(azi));
end

% Remove NaNs
is_nans = isnan(resp);
resp(is_nans) = [];
azi(is_nans) = [];
ele(is_nans) = [];

% initiates each vector's projection onto x,y,z axis respectively
xs=0;     
ys=0;     
zs=0;
  
[azi,sorted] = sort(azi);
resp = resp(sorted);

% Deal with unevenly distributed azimuths.   HH20130422
azi_prime = linspace(min(azi),max(azi),length(azi));  % Not changed if azi is evenly distributed itself. HH20130826
resp = interp1(azi,resp,azi_prime);
azi = azi_prime;


for i=1: length(azi)       
    [xn,yn,zn] = sph2cart(azi(i)*pi/180,ele(i)*pi/180, resp(i) );
    xs = xs + xn;
    ys = ys + yn;
    zs = zs + zn;    
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
vector_amp = vectorsum_amp ;


return;