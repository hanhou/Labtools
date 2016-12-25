function S=EMModeTestCirc(si,fc,r0,m0)
% excess mass mode test
% modified for circular data
%
% si: sampling interval
% fc: function values; circular data
% r0: amplitude threshold
% m0: mass threshold
%
% return
% S: the strength of modality beyond the first (k-1) modes

DEBUG=0;

if DEBUG
    figure(4);
    clf;
    plot([si*(0:length(fc)-1)],fc);
    axis([-si si*length(fc) 0 max(fc)]);
    hold on;
    plot([-si si*length(fc)],[r0 r0],'k');
    pause
end

% total number of points
nfc=length(fc);
Em0=sum(fc)*si; %total mass

%copy the end to the first and the first to the end because fc is circular
fca=zeros(1,nfc+2);
fca(1,1)=fc(1,nfc);
fca(1,2:nfc+1)=fc(1,:);
fca(1,nfc+2)=fc(1,1);

% gradient
gfca=diff(fca); %note: length(gfca)=length(fca)-1
gfca=sign(gfca); % +1: rising, 0: constant, -1: falling

% second gradient; note: length(ggfca)=length(fca)-2
ggfca=diff(gfca); % +2: local minimum, -2: local maximum, +1,-1: edges of flat parts

% add indices
ggfca=cat(1,ggfca,[1:nfc]);
% squeeze out non local min/max points (i.e. zeros)
[i j]=find(ggfca(1,:)~=0);
if ~isempty(j)
    % there is at least one mode
    
    localMinMax=ggfca(:,j);
    % if -1 appears twice in a row, it's a peak with a plateau -> -2
    % if +1 appears twice in a row, it's a trough with a plateau -> +2
    % if -1 is followed by +1 or vice versa, it's a slop with a plateau -> 0
    % find the first peak(-2) or trough(+2)
    [i k]=find(abs(localMinMax(1,:))==2);
    firstMinMax=k(1);
    % starting from the first peak/trough go through the data for one cycle
    for i=0:length(j)-1
        k1=mod(firstMinMax-1+i,length(j))+1;
        k2=mod(firstMinMax+i,length(j))+1;
        if abs(localMinMax(1,k1))==1
            localMinMax(1,k1)=localMinMax(1,k1)+localMinMax(1,k2);
            if localMinMax(1,k1)>0
                localMinMax(1,k2)=localMinMax(1,k1);
            else
                localMinMax(1,k2)=0;
            end
        end
    end
    % eliminate 0
    [i j]=find(localMinMax(1,:)~=0);
    localMinMax=localMinMax(:,j);

    % make lists of minima on the left side of modes, maxima, minima on the right side of modes
    [i j]=find(localMinMax(1,:)<0);
    nmodes=length(j); % number of modes
    modeList=zeros(3,nmodes); % list of modes
    % modeList(1,:): local min on the left side of the mode
    % modeList(2,:): the mode location
    % modeList(3,:): local min on the right side of the mode
    modeList(2,:)=localMinMax(2,j);
    i=length(localMinMax);
    k=mod(j-1+i-1,i)+1;
    modeList(1,:)=localMinMax(2,k);
    k=mod(j,i)+1;
    modeList(3,:)=localMinMax(2,k);
    
    %eliminate modes below r0
    [i j]=find(fc(modeList(2,:))>r0);
    modeList=modeList(:,j);
    nmodes=length(j);
    
    if nmodes>0
        % determine a lamda (mass level) for each mode
        lamda=max(max(fc(modeList(1,:)),fc(modeList(3,:))),r0);

if DEBUG
        for i=1:nmodes
            plot((modeList(2,i)-1)*si,fc(1,modeList(2,i)),'o');
        end
        modeList(:,1:nmodes)
        fc(1,modeList(1,1:nmodes))
        fc(1,modeList(2,1:nmodes))
        fc(1,modeList(3,1:nmodes))
        pause
end

        % compute excess mass
        Em=zeros(1,nmodes);
        a=zeros(1,nmodes);
        b=zeros(1,nmodes);
        for i=1:nmodes
            if modeList(1,i)<modeList(2,i)
                flag=fc(modeList(1,i):modeList(2,i))>=lamda(i);
                [c j]=max(flag);
                a(i)=modeList(1,i)+j-1; % the starting point for integration
            else
                if fc(1)<=lamda(i)
                    flag=fc(1:modeList(2,i))>=lamda(i);
                    [c j]=max(flag);
                    a(i)=j; % the starting point for integration
                else
                    flag=fc(modeList(1,i):nfc)>=lamda(i);
                    [c j]=max(flag);
                    a(i)=modeList(1,i)+j-1; % the starting point for integration
                end
            end
            if modeList(2,i)<modeList(3,i)
                flag=fc(modeList(3,i):-1:modeList(2,i))>=lamda(i);
                [c j]=max(flag);
                b(i)=modeList(3,i)-j+1; % the ending point for integration
            else
                if fc(nfc)<=lamda(i)
                    flag=fc(nfc:-1:modeList(3,i))>=lamda(i);
                    [c j]=max(flag);
                    b(i)=nfc-j+1;
                else
                    flag=fc(modeList(3,i):-1:1)>=lamda(i);
                    [c j]=max(flag);
                    b(i)=modeList(3,i)-j+1;
                end
            end
        
            if a(i)<b(i)
                for j=a(i):b(i)
                    Em(i)=Em(i)+fc(j)-lamda(i); % integrate the area under the function above lamda
if DEBUG
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
end
                end
            else
                for j=a(i):nfc
                    Em(i)=Em(i)+fc(j)-lamda(i); % integrate the area under the function above lamda
if DEBUG
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
end
                end
                for j=1:b(i)
                    Em(i)=Em(i)+fc(j)-lamda(i); % integrate the area under the function above lamda
if DEBUG
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
end
                end
            end
        end
        Em=Em*si; % multiply sampling interval

if DEBUG
        pause
end

        % eliminate modes
        while min(Em)<=m0*Em0 %*Em0 added so that m0 indicates % of total mass
            [c i]=min(Em);

if DEBUG
            if a(i)<b(i)
                for j=a(i):b(i)
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-k');
                end
            else
                for j=a(i):nfc
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-k');
                end
                for j=1:b(i)
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-k');
                end
            end        
end

            if lamda(i)>r0
                nm=length(Em);
                % the mode is not isolated, so it needs to be combined with the neighbor
                if fc(modeList(1,i))==lamda(i)
                    % combine the mode with the one on the left
                    k=mod(i-1+nm-1,nm)+1;
                    modeList(3,k)=modeList(3,i);
                else
                    % combine the mode with the one on the right
                    k=mod(i,nm)+1;
                    modeList(1,k)=modeList(1,i);
                end

                % re-compute excess mass
                % update lamda
                lamda(k)=max(max(fc(modeList(1,k)),fc(modeList(3,k))),r0);
                if modeList(1,k)<modeList(2,k)
                    flag=fc(modeList(1,k):modeList(2,k))>=lamda(k);
                    [c j]=max(flag);
                    a(k)=modeList(1,k)+j-1; % the starting point for integration
                else
                    if fc(1)<=lamda(k)
                        flag=fc(1:modeList(2,k))>=lamda(k);
                        [c j]=max(flag);
                        a(k)=j; % the starting point for integration
                    else
                        flag=fc(modeList(1,k):nfc)>=lamda(k);
                        [c j]=max(flag);
                        a(k)=modeList(1,k)+j-1; % the starting point for integration
                    end
                end
                if modeList(2,k)<modeList(3,k)
                    flag=fc(modeList(3,k):-1:modeList(2,k))>=lamda(k);
                    [c j]=max(flag);
                    b(k)=modeList(3,k)-j+1; % the ending point for integration
                else
                    if fc(nfc)<=lamda(k)
                        flag=fc(nfc:-1:modeList(3,k))>=lamda(k);
                        [c j]=max(flag);
                        b(k)=nfc-j+1;
                    else
                        flag=fc(modeList(3,k):-1:1)>=lamda(k);
                        [c j]=max(flag);
                        b(k)=modeList(3,k)-j+1;
                    end
                end
            
                % update Em
                if a(k)<b(k)
                    Em(k)=0.0;
if DEBUG
                    plot([a(k) b(k)]*si,[lamda(k) lamda(k)],'k');
                    pause
end
                    for j=a(k):b(k)
                        Em(k)=Em(k)+max(fc(j)-lamda(k),0); % integrate the area under the function above lamda
if DEBUG
                        plot([j-1 j-1]*si,[lamda(k) fc(j)],'-r');
end
                    end
                else
if DEBUG
                    plot([a(k) nfc]*si,[lamda(k) lamda(k)],'k');
                    plot([1 b(k)]*si,[lamda(k) lamda(k)],'k');
                    pause
end
                    for j=a(k):nfc
                        Em(k)=Em(k)+max(fc(j)-lamda(k),0); % integrate the area under the function above lamda
if DEBUG
                        plot([j-1 j-1]*si,[lamda(k) fc(j)],'-r');
end
                    end
                    for j=1:b(k)
                        Em(k)=Em(k)+max(fc(j)-lamda(k),0); % integrate the area under the function above lamda
if DEBUG
                        plot([j-1 j-1]*si,[lamda(k) fc(j)],'-r');
end
                    end
                end
                Em(k)=Em(k)*si; % multiply sampling interval
            end

            flag=[1:length(Em)]~=i;
            [i j]=find(flag==1);
            if ~isempty(j)
                Em=Em(j);
                lamda=lamda(j);
                modeList=modeList(:,j);
                a=a(j);
                b=b(j);
        
if DEBUG
                pause
end

            else
                Em=[];
            end
    
        end

        % plot the modes that survived
if DEBUG
        for i=1:length(Em)
            if a(i)<b(i)
                for j=a(i):b(i)
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
                end
            else
                for j=a(i):nfc
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
                end
                for j=1:b(i)
                    plot([j-1 j-1]*si,[lamda(i) fc(j)],'-r');
                end
            end        
        end
end

    else
        Em=[];
    end
    
    % now compute an S metric, the total excess mass beyond k-mode
    if ~isempty(Em)
        % sort excess masses
        sEm=sort(Em);
        % cumulative sum
        cEm=cumsum(sEm);
        % reverse the order
        S=cEm(length(cEm):-1:1);
        % the first entry is the total sum of excess masses,
        % the second entry is the total sum of excess masses beyond one-mode,
        % the third entry is the total sum of excess masses beyond two-modes,
        % etc.
    else
        S=0;
    end

else
    S=0;
end
