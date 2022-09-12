function [f,P,Rot]=MySineSpec(data,samplefreq,K);
%   [f,P]=MySineSpec(data,samplefreq,K) takes an input velocity 
%   vector 'data' (complex), sample frequency in samples per day, K 
%   multitapers (generally 3, and returns a sine multitaper 
%   spectrum P with frequency vector f.   
%   
%   when plotting, loglog(f,P,-f,P) to see both sides.

%   See also MySineTaper (real)  MW 8/03
%   Modified by aleboyer@ucsd.edu 10/09/2021
%  

indg=find(~isnan(data));
data=data(indg);
N=length(data);

if rem(N,2)==0
    data=data(1:end-1);
    N=length(data);
end
u=real(data);
v=imag(data);

u=detrend(u);
v=detrend(v);

rex=complex(u,v);
A=ones(N,K);
for k=1:K
	A(:,k)=((2/(N+1))^(1/2)*sin(((k)*pi*(1:N))/(N+1)))';
	A(:,k)=A(:,k).*rex;  %this gives us h*X in the formula.  
end
fta=fft(A);
msftA=(1/samplefreq)*(abs(fta)).^2;
sumk=sum(msftA,2)./K;  

dz=1/samplefreq;
dk=(1/N)/dz;
%This is not quite right if N is odd...  fix.  mha 6/25/99
if rem(N,2)==0
k=-N/2*dk:dk:N/2*dk-dk;
else
	kp=dk:dk:dk*floor(N/2);
	k=[fliplr(-kp) 0 kp];
end

ip=find(k>0);
im=find(k<0);

f=k;
Pa=fftshift(sumk');
CCW=Pa(ip);
CW=fliplr(Pa(im));
Rot=CW./CCW;
P=0.*Pa;
P(ip)=CCW;
P(im)=fliplr(CW);
