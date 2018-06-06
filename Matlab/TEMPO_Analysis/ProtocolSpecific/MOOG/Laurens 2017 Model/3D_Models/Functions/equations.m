% basic equations
% LBY 20170328

function equations(switch_codes)

switch switch_codes
    case 1 % velocity
        eval('exp((-(t-mu).^2)./(2*sig.^2))');
    case 2 % acceleration
        eval('-(t-mu)./sig.^2.*exp((-(t-mu).^2)./(2*sig.^2))');
    case 3 % jerk
        eval('((t-mu).^2-sig.^2)./sig.^4.*exp((-(t-mu).^2)./(2*sig.^2))');
    case 4 % position
        eval('cumsum(exp((-(t-mu).^2)./(2*sig.^2)))');
end




end