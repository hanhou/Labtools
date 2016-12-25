function Distance=DistanceEvaluate(tempmatrix)
%this function mainly concerns about the sum distance .....

%tempmatrix=z_gauss(2:end-1,2:end-1);
k=0;

for m=1:size(tempmatrix,1)    
    for n=1:size(tempmatrix,2)       
        if tempmatrix(m,n)<max(max(tempmatrix))
            [a,b]=find(tempmatrix==max(max(tempmatrix)));
            tempmatrix(a,b)=-1;
            k=k+1;
            x(k)=a;y(k)=b;            
        else  
        end
    end
end

%Distance=sum(sqrt(diff(x(1:0.5*size(tempmatrix,1)*size(tempmatrix,2))).^2+diff(y(1:0.5*size(tempmatrix,1)*size(tempmatrix,2))).^2));
%Distance=sum(sqrt(diff(x(1:3)).^2+diff(y(1:3)).^2));
Distance=sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2)+sqrt((x(1)-x(3)).^2+(y(1)-y(3)).^2);
