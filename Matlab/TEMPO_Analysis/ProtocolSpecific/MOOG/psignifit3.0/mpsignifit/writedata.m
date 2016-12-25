function writedata ( data, fname='__dummydata.txt' )
% write data to be readible from the command line interface
%
% This file is part of psignifit3 for matlab. (c) 2010 by Ingo Fründ

% Check data format
if size ( data )(2) != 3
    error ( 'data should have three columns' )
end

f = fopen ( fname, 'w' );
for k=1:size ( data )(1)
    fprintf ( f, '%f %d %d\n', data(k,1), data(k,2), data(k,3) );
end
fclose(f);
