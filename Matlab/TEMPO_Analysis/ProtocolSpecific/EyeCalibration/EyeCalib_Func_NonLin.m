function err = EyeCalib_Func_NonLin(a)
% computes the error between the observed and desired eye positions

global Eye_Data

iso_raw = Eye_Data(:,1);
ortho_raw = Eye_Data(:,2);
desired = Eye_Data(:,3);

new = a(1) + a(2).*iso_raw + a(3).*(iso_raw.^2) + a(4).*ortho_raw + a(5).*(ortho_raw.^2);

% for moment, just return the sum squared error
err = norm(new - desired);

