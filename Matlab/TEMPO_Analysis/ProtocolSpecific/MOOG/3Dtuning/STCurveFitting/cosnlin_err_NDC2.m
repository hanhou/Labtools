function err = cosnlin_err_NDC2(vect)
	global rawdata xdata tdata DC2factor
    passdata = [xdata;tdata];    
%     [y,newvect] = funccosnlin(vect,passdata);
    y = funccosnlin_NDC2(vect,passdata);
    y(y<0) = 0;
    err = sqrt(sum( sum(( y - rawdata ) .^2) ));
%     err = sum( sum( sqrt(y) - sqrt(rawdata) ) .^2 );    
return;
 
%     if (DC2factor == 1)
%         y = vonmisesfunc(vect, passdata);
%     elseif (DC2factor == 2)
%         y = vonmisesfuncDC2(vect, passdata);
%     elseif (DC2factor == 3)
%         y = vonmisesfuncDC4(vect, passdata);
% %     elseif (DC2factor == 4)
% %         y = vonmisesfunc2weights(vect, passdata);
%     end
        %     y = y';



