% Smooth kernels for PSTH
% LBY 20170606

function respon = PSTH_smooth( nBins, PSTH_onT, timeWin, timeStep, PSTH_data, kernel_inx, sig)

respon = zeros(nBins,size(PSTH_data,2));

for nn = 1:nBins
    
switch kernel_inx
    case 1
        respon(nn,:) = sum(PSTH_data(PSTH_onT-timeWin/2+timeStep*(nn-1):PSTH_onT+timeWin/2+timeStep*(nn-1),:),1)/(timeWin/1000);
    case 2
        t = -timeWin/2:timeWin/2;
        gau_kernel = exp((-t.^2)./(2*sig.^2));
        gau_kernel = repmat(gau_kernel,size(PSTH_data,2),1);
        gau_kernel = permute(gau_kernel, [2 1]);
        respon(nn,:) = sum(PSTH_data(PSTH_onT-timeWin/2+timeStep*(nn-1):PSTH_onT+timeWin/2+timeStep*(nn-1),:).*gau_kernel,1)/(timeWin/1000);
    case 3
        nan;
end

end

        
        
        
        
% t = -10:10;
% mu = 0;
% sig = 10;
% y = exp((-(t-mu).^2)./(2*sig.^2));
% plot(t,y);

end