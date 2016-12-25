Pd=2000; 
Fd=1; 
Fs=4*Fd; 
R=0.5; 
Delay=5; 
No=1; 
M=4; 
x1=randint(Pd,1,M);  
x2=randint(Pd,1,M); 
x3=randint(Pd,1,M); 
y1=modmap(x1,Fd,Fs,'qask',M); 
y2=modmap(x2,Fd,Fs,'qask',M); 
y3=modmap(x3,Fd,Fs,'qask',M); 
[rcv_a1,ti]=rcosflt(y1,Fd,Fs,'fir/sqrt/Fs',R,Delay) 
[rcv_a2,ti]=rcosflt(y2,Fd,Fs,'fir/sqrt/Fs',R,Delay) 
[rcv_a3,ti]=rcosflt(y3,Fd,Fs,'fir/sqrt/Fs',R,Delay) 
s1=amodce(rcv_a1,10,'qam');  
s2=amodce(rcv_a2,10,'qam'); 
s3=amodce(rcv_a3,10,'qam'); 
 
save sig3   s1 s2 s3 
 
clear 
i=sqrt(-1); 
j=i; 
m=8;  
p=3;    
angle1=10;   
angle2=20; 
angle3=30;  
th=[angle1;angle2;angle3]; 
nn=4096;  
SN1=5;   
SN2=5; 
SN3=5; 
sn=[SN1;SN2;SN3]; 
degrad=pi/180; 
 
 
load sig3 
tt=1:nn; 
S=[s1(tt).';s2(tt).';s3(tt).']; 
nr=randn(m,nn); 
ni=randn(m,nn); 
U=nr+i*ni;   
Ps=S*S'/nn; 
ps=diag(Ps); 
refp=2*10.^(sn/10); 
tmp=sqrt(refp./ps); 
S2=diag(tmp)*S; 
 
tmp=-i*pi*sin(th'*degrad); 
tmp2=[0:m-1]'; 
a2=tmp2*tmp; 
A=exp(a2); 
X=A*S2+U;    
Rxx=X*X'/nn; 
N=m; 
%   numR=size(Rxx,1); 
%   R1=Rxx(1:numR-1,1:numR-1); 
  [V1 D1]=eig(Rxx); 
% D1=sort(D1);  
% D1=D1(end:-1:1);  
   
   
  for n=1:N 
        fenzi=0; 
        fenmu=1; 
        for i=n+1:N 
        fenzi=fenzi+D1(N+1-i:N+1-i,N+1-i:N+1-i); 
        fenmu=fenmu*D1(N+1-i:N+1-i,N+1-i:N+1-i); 
        end 
        fenzi=fenzi/(N-n); 
%         fenmu= nthroot(fenmu,N-n); 
        fenmu=fenmu.^ (1/(N-n)); 
        det(n)=fenzi/fenmu; 
        AICC(n)=2*nn*(N-n)*log(det(n))+2*n*(2*N-n); 
%         fenzi=0; 
%         fenmu=1; 
  end 
     
   
   
   
   
   
%   VV=[V1 zeros(numR-1,1); zeros(1,numR)]; 
%    
%   VV(numR,numR)=1; 
%   R2=VV'*Rxx*VV; 
%   pp=R2(1:numR-1,numR); 
% %   pp=abs(pp); pp=sort(pp); pp=pp(end:-1:1); 
%    
%     for k=1:numR-1 
%         GDE(k)=pp(k)-(1/(numR-1))*mean(pp);  
%         if GDE(k)<=0 & GDE(k-1)>0 
%             signum=k-1;  
%         end 
%     end 