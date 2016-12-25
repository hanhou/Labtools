% 1D azimuth tuning curve
aa = dlmread('globaltuning.TXT');
aa_size = size(aa);

for i = 1 : aa_size(1)
    if (aa(i,1)>270 | aa(i,1)<90) & (aa(i,2)>270 | aa(i,2)<90) % right side of the reference
        sign(i) = 0; % congruent
    elseif (aa(i,1)>270 | aa(i,1)<90) & (aa(i,2)>90 & aa(i,2)<270)
        sign(i) = 180; % opposite
        
    elseif (aa(i,1)>90 & aa(i,1)<270) & (aa(i,2)>90 & aa(i,2)<270) % Left side of the reference
        sign(i) = 0; % congruent
    elseif (aa(i,1)>90 & aa(i,1)<270) & (aa(i,2)<90 | aa(i,2)>270)
        sign(i) = 180; % opposite
    end
end

dlmwrite('globaltuningout.txt',sign');