l = length(FILE);
lpt=length(PATH);

poslpt=lpt;

for i=lpt-1:-1:1
    if (PATH(i)=='\')
        poslpt=i;
        break;
    end
end

if (FILE(l-3:l) == '.htb')	% .htb extension already there
    myfilename = [PATH(1:poslpt) 'Plexon\sorted\' FILE(1:l-3) 'bplx'];
else
    myfilename = [PATH(1:poslpt) 'Plexon\sorted\' FILE '.bplx'];   %the HTB data file
end

if (exist(myfilename,'file'))
    fid=fopen(myfilename,'r');
    fseek(fid,-4*3,'eof');
    mysize=fread(fid,3,'int');
    
    mydata=zeros(mysize(1),mysize(2),mysize(3));
    
    fseek(fid,0,'bof');
   
    for i=1:mysize(1)
        for k=1:mysize(3)
          for j=1:mysize(2)
                mydata(i,j,k)=fread(fid,1,'uchar');
            end
        end
    end
    fclose(fid);
    good_data.spike_data=[good_data.spike_data;mydata];
end

if (FILE(l-3:l) == '.htb')	% .htb extension already there
    myfilename = [PATH(1:poslpt) 'Plexon\sorted\' FILE(1:l-4) '_AD.bplx'];
else
    myfilename = [PATH(1:poslpt) 'Plexon\sorted\' FILE '_AD.bplx'];   %the HTB data file
end

if (exist(myfilename,'file'))
    fid=fopen(myfilename,'r');
    fseek(fid,-4*3,'eof');
    mysize=fread(fid,3,'int');
    
    good_data.plfp_data=zeros(mysize(1),mysize(2),mysize(3));
    
    fseek(fid,0,'bof');
   
    for i=1:mysize(1)
        for k=1:mysize(3)
          for j=1:mysize(2)
                good_data.plfp_data(i,j,k)=fread(fid,1,'short');
            end
        end
    end
    fclose(fid);
end