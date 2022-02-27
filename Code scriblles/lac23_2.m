clc;clear all;close all;

n2=2.8;
n1=1.5;
np=[n1,n2];
d1=160e-9;
d2=80e-9;
% d=50e-9;
lmd=480e-9;
k=2*pi/lmd;
c=3e8;
th1=0;
th2=asind(n1*sind(th1)/n2);
w=k*c; 
a=d1+d2;
k1x=w*n1*cosd(th1)/c;
k2x=w*n2*cosd(th2)/c;

D1=[1,1;n1*cosd(th1),-n1*cosd(th1)];
D2=[1,1;n2*cosd(th2),-n2*cosd(th2)];
d12=[(n2*cosd(th2)+n1*cosd(th1))/2*n2*cosd(th2), (n2*cosd(th2)-n1*cosd(th1))/2*n2*cosd(th2);(n2*cosd(th2)-n1*cosd(th1))/2*n2*cosd(th2),(n2*cosd(th2)+n1*cosd(th1))/2*n2*cosd(th2)];
d21=[(n1*cosd(th1)+n2*cosd(th2))/2*n1*cosd(th1), (n1*cosd(th1)-n2*cosd(th2))/2*n1*cosd(th1);(n1*cosd(th1)-n2*cosd(th2))/2*n1*cosd(th1),(n1*cosd(th1)+n2*cosd(th2))/2*n1*cosd(th1)];
p2=[exp(1i*k2x*d2),0;0,exp(-1i*k2x*d2)];
p1=[exp(1i*k1x*d1),0;0,exp(-1i*k1x*d1)];

m11= exp(-1i*k1x*d1)*(cos(k2x*d2)-1i/2*(k1x/k2x+k2x/k1x)*sin(k2x*d2));
m12= exp(1i*k1x*d1)*(1i/2*(k1x/k2x-k2x/k1x)*sin(k2x*d2));
m21= -exp(-1i*k1x*d1)*(1i/2*(k1x/k2x-k2x/k1x)*sin(k2x*d2));
m22= exp(1i*k1x*d1)*(cos(k2x*d2)+1i/2*(k1x/k2x+k2x/k1x)*sin(k2x*d2));
M1=D1^(-1)*D2*p2*D2^(-1)*D1*p1;
K=1/a*acos((M1(1)+M1(4))/2);

e10=M1(3);
e110=exp(1i*K*a)-M1(1);

% E=(e10*exp(1i*k1x*(x-l*a))+e110*exp(-1i*k1x*(x-l*a)))*exp(-1i*K*l*a);



M2=D2^(-1)*D1*p1*D1^(-1)*D2*p2;
K2=1/a*acos((M2(1)+M2(4))/2);
% e20=M2(3);
% e220=exp(1i*K2*a)-M2(1);
E20=D2^(-1)*D1*p1*[e10;e110];
e20=E20(1);
e220=E20(2);

% e1=1000;
% e11=0;
% %la-d1 to la
% E1=e1*exp(1i*k1x*(x-l*a))+e11*exp(-1i*k1x*(x-l*a));
% %(l-1)a to (l-1)a+d2
% E2=e2*exp(1i*k2x*(x-l*a+d1))+e22*exp(-1i*k1x*(x-l*a+d1));
% 
% out2=d12*p1*[e1;e11];
% out1=d21*p2*[e2;e22];
N=5;
gap=1e-9;
x=N*a;
X=[];
EE=[];
XX=[];
for ii=N:-1:1
    l=ii;
    
    for jj=a/gap:-1:1
        if jj==a/gap
           XX=[XX,x];
        end
        
        if jj>=d2/gap
            
            E1=(e10*exp(1i*k1x*(x-l*a))+e110*exp(-1i*k1x*(x-l*a)))*exp(-1i*K*l*a);
            x=x-gap;
            EE=[EE,abs(E1)];
            X=[X,x];
            
        else
            
            E2=(e20*exp(1i*k2x*(x-(l*a-d1)))+e220*exp(-1i*k2x*(x-(l*a-d1))))*exp(-1i*K2*(l*a-d1));
            x=x-gap;
            EE=[EE,abs(E2)];
            X=[X,x];
            
            
        end
        
    end
       
end


for iii=1:length(XX)
%     figure(1);hold on; xline(XX(iii));
end

figure(1);hold on; plot(X,EE,'r');





