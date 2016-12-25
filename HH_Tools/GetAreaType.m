function type = GetAreaType(inputStr)
global GMTypes;
for i = 1:size(GMTypes,1)
    if sum(findstr(inputStr,GMTypes{i,1})) 
        if ~sum(findstr(inputStr,'?'))
            type = -i;   % Sure
        else   % NOT Sure
            type = -100 - i;
        end
        return
    end
end

% If not find
type = NaN;
return