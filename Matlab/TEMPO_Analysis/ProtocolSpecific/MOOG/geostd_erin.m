% calculate geometric standard deviation, x is a vector 


function y = geostd(x);
x = [2.088
2.574
2.281
2.708
2.689
2.534
2.82
2.199
2.217
2.17
2.818
0.092
2.618
1.996
2.427
1.893
2.166
3.179
1.893
1.729
2.26
2.416
2.505
1.773
];
s = 0;
for i=1:length(x)  
    ss = (log(x(i)) - log(geomean(x)))^2;  
    s = s + ss;
end

ee = sqrt(s/length(x));
y = exp(ee); %geostd




