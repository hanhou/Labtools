% position function -> the integration of Gaussian function
% LBY 20170326

function respon = pos_func(a,t)

mu = a(1);
sig = a(2);

respon = cumsum(exp((-(t-mu).^2)./(2*sig.^2)));
respon = respon./(max(respon)-min(respon));


end