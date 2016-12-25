userName=getenv('UserName');

if (strcmp('Jian',getenv('UserName')) || strcmp('Tatyana Yakusheva',getenv('UserName')))

    if (exist('C:\Documents and Settings\Tatyana Yakusheva','dir')==0)
        mkdir('C:\Documents and Settings\Tatyana Yakusheva');
    end
    if (exist('C:\Documents and Settings\Tatyana Yakusheva\temp','dir')==0)
        mkdir('C:\Documents and Settings\Tatyana Yakusheva\temp');
    end
    ofileP='"C:\Documents and Settings\Tatyana Yakusheva\temp"';

    tmplen=max(size(FILE));

    FILEp=FILE(1:tmplen-4);

    lpt=length(PATH);

    if PATH(lpt)~='\'
        PATH=[PATH '\'];
    end

    lpt=length(PATH);

    for i=lpt-2:-1:1
        if (PATH(i)=='\')
            poslpt=i;
            break;
        end
    end

    %%%Set path name for the location where smr files are stored
    %%% If the SMRs are located at "CED", uncomment the following line
    %     PATH1=[PATH(1:poslpt) 'CED\'];

    %%% If the SMRs are located at "smr_sorted", uncomment the following line
    PATH1='Z:\Data\MOOG\Priam\Analysis\Tanya\smr_sorted\';

    %%% If the spike channel name is "Spike1", uncomment the following line
%         cmd=sprintf('recoverspk -p %s -i %s -o %s -l %s',PATH1,FILEp,ofileP,PATH);

    %     cmd=sprintf('recoverspk -p %s -i %s -o %s -l %s -s Spike1',PATH1,FILEp,ofileP,PATH);
    %%% Otherwise, uncomment the following line and then replace SpikeName with an actual spike channel name
    cmd=sprintf('recoverspk -p %s -i %s -o %s -l %s -s Spikes -e marker',PATH1,FILEp,ofileP,PATH);
    



    dos(cmd);

    l = length(FILE);

    lpt2=length(ofileP);
    myfilename = [ofileP(2:lpt2-1) '\' FILE(1:l-3) 'spk'];

    if (exist(myfilename,'file'))

        clear data.spike_data;

        fid=fopen(myfilename,'r');

        fseek(fid,-4*3,'eof');
        mysize=fread(fid,3,'int');

        data.spike_data=zeros(2,mysize(2),mysize(3));

        fseek(fid,0,'bof');

        for i=1:mysize(1)
            for k=1:mysize(3)
                for j=1:mysize(2)
                    data.spike_data(i,j,k)=fread(fid,1,'uchar');
                end
            end
        end

        fclose(fid);

        for k=1:mysize(3)
            for j=1:1000
                data.spike_data(1,j,k)=data.spike_data(1,j+4000,k);
            end
        end

        for k=1:mysize(3)
            for j=1017:17:3000
                data.spike_data(2,j,k)=1;
            end
        end

        %markers=[4,5,12,13,15];
        %idx=zeros(1,5);
        %for k=1:mysize(3)
        %    for j=1:5
        %       idx(j)=find(data.event_data(1,:,k)==markers(j));
        %    end
        %
        %    sidx=idx(1);
        %    for j=2:5
        %        data.event_data(1,idx(j),k)=0;
        %        nidx=sidx+1015;
        %        if nidx>5000
        %            nidx=5000;
        %        end
        %        data.event_data(1,nidx,k)=markers(j);
        %    end
        %end
    end

    myfilename = [ofileP(2:lpt2-1) '\' FILE(1:l-3) 'evn'];

    if (exist(myfilename,'file'))

        clear data.event_data;

        fid=fopen(myfilename,'r');

        fseek(fid,-4*3,'eof');
        mysize=fread(fid,3,'int');

        data.event_data=zeros(mysize(1),mysize(2),mysize(3));

        fseek(fid,0,'bof');

        for i=1:mysize(1)
            for k=1:mysize(3)
                for j=1:mysize(2)
                    data.event_data(i,j,k)=fread(fid,1,'uchar');
                end
            end
        end

        fclose(fid);

    end

end

