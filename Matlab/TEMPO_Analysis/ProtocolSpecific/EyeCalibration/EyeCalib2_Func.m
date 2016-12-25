function err = EyeCalib2_Func(a)
% computes the error between the observed and desired eye positions

global Eye_Data

iso_raw = Eye_Data(:,1);
ortho_raw = Eye_Data(:,2);
desired = Eye_Data(:,3);

new = a(1) + a(2).*iso_raw + a(3).*ortho_raw + a(4).*iso_raw.*ortho_raw;

% for moment, just return the sum squared error
err = norm(new - desired);

