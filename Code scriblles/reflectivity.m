close all;clear all;
% % PERIODICLAYERED MEDIA
% 
% ######### BRAGG REFLECTION ########
% same concept N layers and only care about output and input not the nth layer
% A  B  ^N         AUNm-UNmm          BUNm
%             =   
% C  D             CUNm             DUNm-UNmm
% 
% UN= sin(N+1)*kz*v/sin(kz*v)
% 
% rN=CUNm/(AUNm-UNmm)
% 
% in the form of moldulus
% 
% mod(rN)=sqrt(mod(C)^2/((mod(C)^2)+(sin(kzv)/sin(Nkzv)^2)))

% Define layers and refective constants
n2=2.8;
n1=1.5;
np=[n1,n2];
a=80e-9;
b=160e-9;

% numbere of layers n
N=10;


th=0;
v=a+b;
% z=0;
lmd=(100:0.01:1000)*1e-9;
% ww=(0:0.01:pi);
for i=1:length(lmd)
    k=2*pi/lmd(i);
    c=3e8;
    % form bragg condition
    w=k*c;
%     w=ww(i);
%     k=w/c;
%     if mod(i,2) == 0
% %         w=k*c/(2*pi*np(2));
%         ky=k*np(2)*sin(th);
%     else
% %         w=k*c/(2*pi*np(1));
    ky1=k*np(1)*sind(th);  
    ky2=k*np(2)*sind(th);  
%     end
    k1z=sqrt((np(1)*w/c)^2-ky1^2);
    k2z=sqrt((np(2)*w/c)^2-ky2^2);
    A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    C= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
    D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    K=1/v*acos((A+D)/2);
    UN= sin((N+1)*K*v)/sin(K*v);
    UNm=sin((N)*K*v)/sin(K*v);
    UNmm=sin((N-1)*K*v)/sin(K*v);
    % coeff. of reflectance
    rN(i)=abs(C*UNm/(A*UNm-UNmm));
    % reflectivity
%     rN(i)=sqrt(abs(C)^2/(abs(C)^2+(sin(K*v)/sin(N*K*v))^2));
    %transmitivity
    tN(i)=1-rN(i);
    lambda(i)=lmd(i); 
end
figure(1);hold on; plot(lambda,rN,'r');
% figure(1);hold on; plot(lambda,tN,'b');




