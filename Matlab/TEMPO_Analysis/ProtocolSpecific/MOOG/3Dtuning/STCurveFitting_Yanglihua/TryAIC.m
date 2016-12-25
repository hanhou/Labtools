clear all
x = [ 53 57 58 63 66 67 68 69 70 72 73 75 76 78 79 81 ]';
y = [ 1 1 1 1 0 0 0 0 2 0 0 1 0 0 0 0 ]';
n = [ 1 1 1 1 1 3 1 1 4 1 1 2 2 1 1 1 ]';

[betas dev stats]= glmfit(x,[y n],'binomial', 'link', 'logit');

fitted = glmval(betas,x,'logit');

%% LL should be -10.1576
logLikelihood = sum(log( binopdf( y, n,fitted)))

%% AIC should be 24.3152
AIC = -2*logLikelihood + 2*numel(betas)

 - gD
 
% Background of AIC
%    AIC is used to compare the quality of nested models
%   AIC requires a bias-adjustment in small sample sizes

 
 
%AICad Function
%Date:Dec.29th, 2008 
%Purpose: Calculate the adjusted 
%               Akaike¡¯s Information Criteria (AIC)
%InPut: AIC 
%          n (the sample size)
%          K (the number of parameters in the model)
%OutPut: AICad
%              delta
%              weight
%Functions Called: None
%Called by: None

 
 
>>AIC=[-123; -241.7; -92.4; -256.9]; 
>>n=12;
>>K=[4;8;2;6];
>>[AICad, delta, w] = CalAIC(AIC, n, K)
AICad =            delta =                   w =

  -83.0000              89.9000                   0.0000
  -97.7000              75.2000                   0.0000
  -80.4000              92.5000                   0.0000
 -172.9000                        0                   1.0000


 function[AICad,delta,w]=CalAIC(AIC,n,K)
  AICad = AIC+2.*K.*(K+1);
  a = min(AICad);
  delta = AICad - a;
  b = exp(-0.5.*delta);
  c = sum(b);
  w = b./c;

  
  
  