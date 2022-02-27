close all;clear all;

n1=1.5;
n2=2.8;
np=[n1,n2];
a=100e-9;
b=150e-9;

lmd=500e-9;
k=2*pi/lmd;
c=3e8;

w=k*c; 

% at angle theta ky=k*sin(th)
th=60;

% numbere of layers n
n=50;


v=a+b;
z=0;

for i=0:n
    ky=k*np(1)*sind(th);
    k1z=sqrt((np(1)*w/c)^2-ky^2);
    k2z=sqrt((np(2)*w/c)^2-ky^2);
    A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
    D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    
    K=1/v*acos((A+D)/2);
    a0=B;
    b0=exp(1i*K*v)-A;
    
    E=abs(((a0*exp(-1i*k1z*(z-i*v))+b0*exp(1i*k1z*(z-i*v)))*exp(1i*K*(z-i*v)))*exp(-1i*K*z));
    
    Z(i+1)=z;
    EE(i+1)=E;    
    z=z+10*1e-9;
end
figure(1);hold on; plot(Z,EE,'r');





