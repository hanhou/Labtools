% qphase2.m
%
% Attempts to compress a set of phase values with a large
% spread (>2pi) into a set close to [-pi,pi]

function [qph]=qphase(ph)

qph=sign(ph).*(abs(ph)-floor(abs(ph)/360)*360);
qph=qph + (qph<-180)*360 - (qph>180)*360;
qph=180/pi*unwrap(pi/180*qph);

return


%qph=sign(ph).*(abs(ph)-floor(abs(ph)/360)*360);
%qph=qph + (qph<-180)*360 - (qph>180)*360;
%qph=180/pi*unwrap(pi/180*qph);
