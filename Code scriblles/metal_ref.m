clear all;clc;
n2=1.5;
n1=2.8;
np=[n1,n2];
a=180e-9;
b=90e-9;
d=50e-9;
% at angle theta ky=k*sin(th)
th=60;

v=a+b;
z=0;
N=10;

lmd=(350:0.01:800)*1e-9;

for i=1:length(lmd)
    outm=m1(lmd(i),th,d);
    k=2*pi/lmd(i);
    c=3e8;
    w=k*c;

    ky1=k*np(1)*sind(th);  
    ky2=k*np(2)*sind(th);  
    k1z=sqrt((np(1)*w/c)^2-ky1^2);
    k2z=sqrt((np(2)*w/c)^2-ky2^2);
    A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
    C= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
    D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
    K=1/v*acos((A+D)/2);
    UN= sin((N+1)*K*v)/sin(K*v);
    UNm=sin((N)*K*v)/sin(K*v);
    UNmm=sin((N-1)*K*v)/sin(K*v);
    bm=[A*UNm-UNmm,B*UNm;C*UNm,D*UNm-UNmm];
    a0b0=outm*bm*[1;0];
    rN(i)=abs(a0b0(2)/a0b0(1));
    %transmitivity
    tN(i)=1-rN(i);
    lambda(i)=lmd(i);  
end
figure(1);hold on; plot(lambda,rN,'k');
% figure(1);hold on; plot(lambda,tN,'b');









