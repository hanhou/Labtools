function p=cum_gaussfit(beta,x)
MU=beta(1);
SIGMA=beta(2);
p = normcdf(x,MU,SIGMA);