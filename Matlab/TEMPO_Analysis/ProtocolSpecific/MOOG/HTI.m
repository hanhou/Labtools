% Calculate Heading Tuning Index, need parameters of raw data matrix,
% spontaneous activity

function r = HTI(resp, spon)

if abs( max(max(resp))-spon) > 1      % avoid divied by 0
     M_norm = (resp - spon) / ( max(max(resp))-spon); %calculate normalized vector first 
else
     M_norm = resp - spon;
end

if size(resp) == [5 8] % must omit for 1D case. CRF
    M_norm(1,2:8)=0;
    M_norm(5,2:8)=0;                                                                 % hard code temperarilly
end

[azi_norm,ele_norm,amp_norm] = vectorsum(M_norm);
Vec_sum_norm = [azi_norm,ele_norm,amp_norm];                                         % calculate amplitude of vectorsum
r = abs(Vec_sum_norm(3)) / sum(sum(abs(M_norm))) ;                                   % this is index of spacial focus 
