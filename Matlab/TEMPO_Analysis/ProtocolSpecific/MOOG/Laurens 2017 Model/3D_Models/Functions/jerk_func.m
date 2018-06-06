% jerk function -> the Diff of Diff of Gaussian function
% LBY 20170326

function respon = jerk_func(a,t)

mu = a(1);
sig = a(2);

respon = ((t-mu).^2-sig.^2)./sig.^4.*exp((-(t-mu).^2)./(2*sig.^2));
respon = respon./(max(respon)-min(respon));

end