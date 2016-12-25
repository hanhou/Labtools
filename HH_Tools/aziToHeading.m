function azi = aziToHeading(x)
azi = mod(270-x,360)-180; % Transferred from azimuth to heading
