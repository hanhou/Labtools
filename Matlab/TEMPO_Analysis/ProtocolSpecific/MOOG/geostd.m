% calculate geometric standard deviation, x is a vector
% 09/2006, GY 

function y = geostd(x);

dim = size(x);

for j = 1:dim(2) % column
    s=0;
    for i = 1:dim(1) % row
        ss = ( log(x(i,j)) - log(geomean(x(:,j))) ) ^2;
        s = s + ss;
    end
    ee = sqrt(s/dim(1));
    y(j) = exp(ee); %geostd
end


%see http://en.wikipedia.org/wiki/Geometric_standard_deviation
