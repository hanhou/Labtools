function [mode, pval]=modalityTestForYong(xh,modelist,Nb)
%function [mode, pval]=modalityTestForYong(xh,modelist,Nb)
% This program conducts an excess mass mode test
% on a circular distribution using an algorithm by Fisher & Marron (2001).
% 
% xh: a vector of input data
% modelist: a vactor that lists the numbers of modes to test for
% Nb: a number of bootstrap simulations
%
% return
% mode: the number of modes
% pval: a vector of p-values for each of the numbers of modes tested
%
% modified for Yong. 9/25/06
% main routine

DEBUG=0;

xhrad=pi*xh/180.0;
xhrad=sort(xhrad);

% Step 1: obtain a k-modal kernel density estimate
% compute a kernel density estimate
n=length(xh); % the number of data
nm=length(modelist); % numbers of modes to examine
densitySampleSize=1000;
densitySampleInterval=2*pi/densitySampleSize;
xk=0:densitySampleInterval:2*pi;
r0=0.0; % amplitude threshold
m0=0.0; % excess mass threshold

pval=zeros(1,nm);
mode=0;
for j=1:nm

	if DEBUG
        figure(1);
        clf;
        orient landscape;
        subplot('Position',[0.2 0.15 0.6 0.05]);
        plot([0 2*pi],[0 0],'k');
        hold on;
        plot(xhrad,zeros(1,n),'x');
        hold off;
        axis([0 2*pi 0 1]);
        axis off;
	end

	kl=0.25; % lower bound
	ku=500; % upper bound
	while (ku-kl)>0.0001*kl %search the minimum window size until the difference between the upper bound and the lower bound is less than 1% of the lower bound
        k=(ku+kl)/2; %kernel window size
        
        % compute a kernel density estimate
        Fk=compKDEvonMises(xk,xhrad,k);
        h=acos(1-1/(2*k)); %this is my method
        if DEBUG
            subplot('Position',[0.2 0.25 0.6 0.3]);
            plot(xk,Fk*densitySampleInterval,'-r');
            axis([0 2*pi 0 max(Fk)*densitySampleInterval]);
            pause
        end
        S=EMModeTestCirc(densitySampleInterval,Fk,r0,m0);
	
        %check if Fk is nm-modal
        if length(S)>modelist(j) & S(modelist(j)+1)~=0
            ku=k; % more than k modes
        else
            kl=k; % no more than k modes
        end;
	end
	k=kl;
	Fk=compKDEvonMises(xk,xhrad,k);
	h=acos(1-1/(2*k)); %this is my method
	if DEBUG
        plot(xk,Fk*densitySampleInterval,'-k');
        axis([0 2*pi 0 max(Fk)*densitySampleInterval]);
        pause
	end
	S=EMModeTestCirc(densitySampleInterval,Fk,r0,m0);
	k0=k;
	
	% plot the kernel density estimate
	
	% Step 2: compute a goodness-of-fit test statistic
	%         for the kernel density estimate
	cumFk=cumsum(Fk);
	cumFk=cumFk./cumFk(end);
	Fh=cumFk(floor(xhrad/densitySampleInterval+0.5)+1);
	
	if DEBUG
        subplot('Position',[0.2 0.6 0.6 0.3]);
        plot(1:n,Fh,'-k');
        axis([0 n 0 1]);
        hold on
        plot([1 n],[0 1],'-r');
        pause    
	end
	
	% compute a Watson's U2 statistic
	U2N=WatsonU2N(Fh);
	
	% Step 3: Repeating steps 1 and 2 for bootstrap samples
	U2bootstrap=zeros(1,Nb);
	for b=1:Nb
        
        % generate random samples based on Fk
        Xb=genRandDiscVar(Fk, n, 0);
        Fb=cumFk(sort(Xb));
        
        if DEBUG & mod(b,100)==1
            plot(1:n,Fb,'-g');
            pause
        end
        
        %compute U2 statistic
        U2bootstrap(b)=WatsonU2N(Fb);
	
	end
	
	% Step 4: Computing a significance probability for the goodness-of-fit
	%         test statistic of the original k-modal kernel density estimate
	
	p=U2bootstrap-U2N>0;
	pval(1,j)=sum(p)/Nb;
	
	if pval(1,j)>0.05
        if mode==0
            mode=modelist(j);
        end
	end
    
end
