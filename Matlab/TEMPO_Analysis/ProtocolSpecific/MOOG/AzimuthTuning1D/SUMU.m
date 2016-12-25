ssu=dlmread('SU.txt'); %reads in singleunit data
msu=dlmread('MU.txt'); %reads in multiunit data
dim=size(ssu); %determines the number of neurons to analyze

for i=1: dim(1) %one loop for each cell 

    %first we find the correlation coefficient and p-value for the vestibular condition
   [rr,pp] = corrcoef(ssu(i, 2:11),msu(i, 2:11)); %uses data from relevant columns in text files
   r_ves(i) =rr(1,2); %correlation coefficient
   p_ves(i)= pp(1,2); %p-value
   
   %then we find the correlation coefficient and p-value for the visual condition
   [rr,pp] = corrcoef(ssu(i, 12:21),msu(i, 12:21)); %uses data from relevant columns in text files
   r_vis(i) =rr(1,2); %correlation coefficient
   p_vis(i)= pp(1,2); %p-value

end  %all done :)