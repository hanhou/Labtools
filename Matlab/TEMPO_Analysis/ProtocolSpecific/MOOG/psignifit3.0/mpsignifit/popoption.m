function [value,args] = popoption ( args )
% pop the first value of args and return args
%
%
% This file is part of psignifit3 for matlab (c) by Ingo Fründ

value = args{1};
if length(args)>1
    args = args(2:end);
else
    args = {};
end
