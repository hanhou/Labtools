function data = XlsCell2Mat (raw)
% Turn Excel cells data to matrices.  HH20140604
% Output: data = {matrix(column 1), matrix(column 2), ...}

% Preallocation
[rows, cols] = size(raw);
data = cell(1,cols);

for col = 1:cols
    
    mat = ones(rows,100) * NaN;  % Temporally preallocation
    maxLen = 0;
    for row = 1:rows
        if ~isnan(raw{row,col})
            temp = str2num(raw{row,col});
        else
            temp = NaN;
        end
        
        if sum(isnan(temp)) || isempty(temp)  % Exceptions
            mat(row,1) = NaN;  % Mark NaNs
        else
            mat(row,1:length(temp))=temp;
            if length(temp)>maxLen ;maxLen = length(temp); end
        end
    end

    % Fill NaNs and clear up
    if maxLen > 0 
        mat(:,maxLen+1:100) = [];
    else % All NaNs
    end
    
    mat(isnan(mat(:,1)),1:maxLen) = NaN;
    data{col} = mat;
end

end