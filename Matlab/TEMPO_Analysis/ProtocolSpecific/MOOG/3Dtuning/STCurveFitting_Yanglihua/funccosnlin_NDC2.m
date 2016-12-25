function z= vonmisesfunc(vect,xtdata)%[z,newvect]= vonmisesfunc(vect,xtdata)
    xdata=xtdata(1:9,:);   %      xdata = xtdata(1:8,:);
     tdata=xtdata(10:18,:);%     tdata = xtdata(9:16,:);
    global model_use;
    if (model_use == 1)                %%% func = A * [F(theta) * G(tau)]
        % vect = [Ro  Amplitude n mu_Azi  mu_T  sigma_t];
        %bound = [0-peaktotrough 0-peak 0-5 0-360 0-2 0-2 0-1 -0.5~0];
        gauss_time = exp((-(tdata-vect(5)).^2) / (2*vect(6).^2)) ;
        gauss_time = gauss_time / max(max(gauss_time)); % Normalize Gaussian
        R2 = [[exp(vect(3) * cos(xdata - vect(4)) )] - 1] / vect(3);         
        space_gauss = (R2 - min(min(R2)))/(max(max(R2)) - min(min(R2)));   
        func = vect(2) * [ space_gauss .* gauss_time];
	    z = vect(1) + func;    
    elseif (model_use == 2)            
        %%% func = A * [w* F(theta) * G(tau) + (1 - w) * F (theta + theta1)* dG{tau}]+DC        
        % vect = [Ro  Amplitude n mu_Azi  mu_T  sigma_t wVel ThetaAcc];
        gauss_time = exp((-(tdata-vect(5)).^2) / (2*vect(6).^2)) ;
        for i = 1:size(gauss_time,1)
            d_gauss_time(i,:) = diff(gauss_time(i,:));    
        end
        d_gauss_time(:,end+1) = d_gauss_time(:,end);
        gauss_time = gauss_time / max(max(gauss_time)); % Normalize Gaussian
        if max(max(d_gauss_time))==0
            d_gauss_time = d_gauss_time / (max(max(d_gauss_time))+0.00001);  % Normalize derivative 
        else
             d_gauss_time = d_gauss_time / max(max(d_gauss_time));  % Normalize derivative
        end
        R2 = [exp(vect(3) * cos((xdata - vect(4)))) - 1] ./ vect(3);      
        space_gauss = (R2 - min(min(R2)))/(max(max(R2)) - min(min(R2)));         
        
        R3 = [exp(vect(3) * cos((xdata - vect(4) - vect(8)))) - 1] ./ vect(3);                
        space_gauss_1 =  (R3 - min(min(R3)))/(max(max(R3)) - min(min(R3)));         
        func = vect(2) * [ (vect(7) * (space_gauss .* gauss_time)) + ((1-vect(7)) * (space_gauss_1 .* d_gauss_time)) ];                    
        z = vect(1) + func;
        
     elseif(model_use == 3)         %%% func = A * [w* F(theta) * G(tau) + (1 - w) * F (theta + theta1)* dG{tau}]
        % vect = [Ro  Amplitude  n  mu_Azi  mu_T  sigma_t DC2 Weight1 Theta1 Weight2 Theta2 Weight3];
        gauss_time = exp((-(tdata-vect(5)).^2) / (2*vect(6).^2)) ;
        for i = 1:size(gauss_time,1)
            d_gauss_time(i,:) = diff(gauss_time(i,:));               
            for j=2:size(gauss_time,2)
                clear dt;dt=tdata(i,1:j);
                clear dVel; dVel=gauss_time(i,1:j);
                Pos_gauss_time(i,j)= trapz(dt,dVel);    
            end              
        end        
        d_gauss_time(:,end+1) = d_gauss_time(:,end);
        gauss_time = gauss_time / max(max(gauss_time)); % Normalize Gaussian
        if max(max(d_gauss_time))==0
            d_gauss_time = d_gauss_time / (max(max(d_gauss_time))+0.00001);  % Normalize derivative 
        else
             d_gauss_time = d_gauss_time / max(max(d_gauss_time));  % Normalize derivative
        end   
        Pos_gauss_time = Pos_gauss_time / max(max(Pos_gauss_time));  % Normalize Integration
        R2 = [exp(vect(3) * cos((xdata - vect(4)))) - 1] ./ vect(3);  
        space_gauss = (R2 - min(min(R2)))/(max(max(R2)) - min(min(R2))); 
        
        R3 = [exp(vect(3) * cos((xdata - vect(4) - vect(8)))) - 1] ./ vect(3);        
        space_gauss_1 = (R3 - min(min(R3)))/(max(max(R3)) - min(min(R3))); 
        
        R4=[exp(vect(3)*cos(xdata-vect(4)-vect(9)))-1]./vect(3);      
        space_gauss_2 = (R4 - min(min(R4)))/(max(max(R4)) - min(min(R4))); 
        
        func = vect(2) * [ vect(7) * (space_gauss .*gauss_time) + (1-vect(7)) * (space_gauss_1 .* d_gauss_time) + vect(10)*(space_gauss_2.*Pos_gauss_time)];         
        z = vect(1) + func;  
    end
     z(z<0) = 0;
return;

% F = x(1) * exp(-2*(1-cos(xdata-x(2)))/(x(3))^2) + x(4);