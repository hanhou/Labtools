function err = EyeCalib3_Func(a)
% computes the error between the observed and desired eye positions

global Eye_Data

horiz_raw = Eye_Data(:,1);
vert_raw = Eye_Data(:,2);
horiz_desired = Eye_Data(:,3);
vert_desired = Eye_Data(:,4);

horiz_new = a(1) + a(2).*horiz_raw + a(3).*vert_raw + a(4).*horiz_raw.*vert_raw;
vert_new = a(5) + a(6).*vert_raw + a(7).*horiz_raw + a(8).*horiz_raw.*vert_raw;

% for moment, just return the sum squared error
%err = norm(horiz_new - horiz_desired) norm(horiz_new - horiz_desired);
err = norm([horiz_new; vert_new] - [horiz_desired; vert_desired]);

