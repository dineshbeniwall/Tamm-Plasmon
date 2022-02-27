close all;clear all;
% % PERIODICLAYERED MEDIA
% 
% n(z) = n2, 0<z<b
%        n1, b<z<v
% 
% a  b    a      b       a     b  a  b  a  b 
% n1 n2   n1     n2      n1    n2 n1 n2 n1 n2
%        anm1   an2     an1
%        bnm1   bn2     bn1
% 
% 
%        
% general solution of wave equation can be given by  E(z)*exp(1i*w*t-ky*y)
% electric field distribution in the same layer can be written as :
% 
% E(y,z)=[anp*exp()-1i*kaz*(z-nv)+bnp*exp(1i*kaz*(z-nv))]*exp(-1i*ky*y)
% 
% with kpz=sqrt((npw/c)^2-ky^2)
% 
% anm       A   B  an
%       =   
% bnm       C   D  bn
% 
% where: TE modes
% 
% A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
% 
% B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
% 
% C= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
% 
% D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
% 
% so In case of n layers
% a0      A    B  ^n   an
%     =   
% b0      C    D       bn
% 
% In case of other way around
% an      D   -B  ^n   a0
%     =
% bn     -C    A       b0
% 
% ########  Bloch waves and band structures ###########
% E = Ek(z)*exp(-1i*kz)*exp(1i*(w*t-ky*y))
% where
% Ek(z)= Ek(z+v)
% 
% form above eqn:
% a0        B
%     =   
% b0     exp(1i*K*v)-A
% 
%here K is bloch wave number
% and 
% 
% an                            B
%     =   exp(-1i*n*kz*v)
% bn                        exp(1i*kz*v)-A
% 
% Kz for the bloch wave function
% 
% k(ky,w)=1/v*acos((A+D)/2)
% 
% # plot bnad gap form above eqn between w vs ky or k
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
n1=1.5;
n2=2.8;
np=[n1,n2];
a=100e-9;
b=150e-9;
%k=2pi/lambda
lmd=500e-9;
k=2*pi/lmd;
c=3e8;
% form bragg condition
w=k*c; %for once

% at angle theta ky=k*sin(th)
th=0;

% numbere of layers n
n=30;



k1z=sqrt((np(1)*w/c)^2-ky^2);
k2z=sqrt((np(2)*w/c)^2-ky^2);


% for here we have taken only two slabs
A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));

B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));

C= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));

D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));

v=a+b;
z=0;

for i=1:n
    K=1/v*acos((A+D)/2);
    a0=B;
    b0=exp(1i*k1z*v)-A;
    outn=[D,-B;-C,A]^i*[a0;b0];
    if mod(i,2) == 0
        ky=k*np(2)*sin(th);
        z=z+b;
        E=abs(outn(1)*exp(-1i*k2z*(z-n*v))+outn(2)*exp(1i*k2z*(z-n*v)));
    else 
        z=z+a;
        E=abs(outn(1)*exp(-1i*k1z*(z-n*v))+outn(2)*exp(1i*k1z*(z-n*v)));
    end
%     figure(1);hold on; plot(z,E,'r.');
    Z(i)=z;
    EE(i)=E;    
end
figure(1);hold on; plot(Z,EE,'r');





