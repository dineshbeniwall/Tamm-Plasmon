close all;clear all;

n1=1.5;
n2=2.8;
np=[n1,n2];
a=100e-9;
b=150e-9;

% numbere of layers n
N=30;

z=0;
th=60;
v=a+b;
% z=0;
lmd=(100:0.1:1000)*1e-9;
% ww=(0:0.01:pi);

for i=1:length(lmd)
    k=2*pi/lmd(i);
    c=3e8;
    % form bragg condition
    w=k*c;
    ky=k*np(1)*sind(th);
    k1z=sqrt((np(1)*w/c)^2-ky^2);
    k2z=sqrt((np(2)*w/c)^2-ky^2);
    A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
    D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    
    K=1/v*acos((A+D)/2);
    kk(i+1)=K;
    kyy(i+1)=ky;
    ww(i+1)=w;
    z=z+v;
end
figure(1);hold on; plot3(kk,kyy,ww,'r');xlabel('ky');ylabel('w');zlabel('KK');
% surf(kyy,w,kk);




