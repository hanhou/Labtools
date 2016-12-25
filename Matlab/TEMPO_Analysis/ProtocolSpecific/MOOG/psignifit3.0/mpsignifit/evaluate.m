function predicted = evaluate(x, inference)
% Evaluate the psychometric function model at positions given by x
%
% x 
%   Values at which PF should be evaluated
%
% inference
%   The inference object

inputData = [x; zeros(1, length(x)); zeros(1, length(x))]'

if inference.gammaislambda
    diag = Diagnostics ( inputData, inference.params_estimate, ...
        'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts, 'gammaislambda' );
else
    diag = Diagnostics ( inputData, inference.params_estimate, ...
        'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts );
end

predicted = diag.prediction(:,2)'