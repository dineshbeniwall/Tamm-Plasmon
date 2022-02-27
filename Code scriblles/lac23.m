clc;clear all;close all;

n2=2.8;
n1=1.5;
np=[n1,n2];
d1=80e-9;
d2=180e-9;
% d=50e-9;
lmd=600e-9;
k=2*pi/lmd;
c=3e8;
th1=0;
th2=asind(n1*sind(th1)/n2);
w=k*c; 
a=d1+d2;
k1x=w*n1*cosd(th1)/c;
k2x=w*n2*cosd(th2)/c;

% m11= exp(-1i*k1x*d1)*(cos(k2x*d2)-1i/2*(k1x/k2x+k2x/k1x)*sin(k2x*d2));
% m12= exp(1i*k1x*d1)*(1i/2*(k1x/k2x-k2x/k1x)*sin(k2x*d2));
% m21= -exp(-1i*k1x*d1)*(1i/2*(k1x/k2x-k2x/k1x)*sin(k2x*d2));
% m22= exp(1i*k1x*d1)*(cos(k2x*d2)+1i/2*(k1x/k2x+k2x/k1x)*sin(k2x*d2));

D1=[1,1;n1*cosd(th1),-n1*cosd(th1)];
D2=[1,1;n2*cosd(th2),-n2*cosd(th2)];
d12=[(n2*cosd(th2)+n1*cosd(th1))/2*n2*cosd(th2), (n2*cosd(th2)-n1*cosd(th1))/2*n2*cosd(th2);(n2*cosd(th2)-n1*cosd(th1))/2*n2*cosd(th2),(n2*cosd(th2)+n1*cosd(th1))/2*n2*cosd(th2)];
d21=[(n1*cosd(th1)+n2*cosd(th2))/2*n1*cosd(th1), (n1*cosd(th1)-n2*cosd(th2))/2*n1*cosd(th1);(n1*cosd(th1)-n2*cosd(th2))/2*n1*cosd(th1),(n1*cosd(th1)+n2*cosd(th2))/2*n1*cosd(th1)];
p2=[exp(1i*k2x*d2),0;0,exp(-1i*k2x*d2)];
p1=[exp(1i*k1x*d1),0;0,exp(-1i*k1x*d1)];

e1=1;
e11=0;
% %la-d1 to la
% E1=e1*exp(1i*k1x*(x-l*a))+e11*exp(-1i*k1x*(x-l*a));
% %(l-1)a to (l-1)a+d2
% E2=e2*exp(1i*k2x*(x-l*a+d1))+e22*exp(-1i*k1x*(x-l*a+d1));
% 
% out2=d12*p1*[e1;e11];
% out1=d21*p2*[e2;e22];

gap=0.1e-9;
N=4;
% x=(N*a:gap:0);
x=N*a;
XX=[];
EE=[];
X=[];
for ii=N:-1:1
    l=ii;
    
    for jj=a/gap:-1:1
        
        if jj==a/gap
            
            if l==N
                XX=[XX,x];
                
            else
                out1=D1^(-1)*D2*p2*[e2;e22];
                e1=out1(1);
                e11=out1(2);
                XX=[XX,x];
                
            end
        elseif jj==d2/gap
            
            out2=D2^(-1)*D1*p1*[e1;e11];
            e2=out2(1);
            e22=out2(2);
            XX=[XX,x];
        end
        
        if jj>d2/gap
            E1=e1*exp(1i*k1x*(x-l*a))+e11*exp(-1i*k1x*(x-l*a));
%             E1=e1*exp(1i*k1x*(x-l*a));
            X=[X,x];
            EE=[EE,abs(E1)];
            x=x-gap;
            
            
        else
            E2=e2*exp(1i*k2x*(x-l*a+d1))+e22*exp(-1i*k2x*(x-l*a+d1));
%             E2=e2*exp(1i*k2x*(x-l*a+d1));
            X=[X,x];
            EE=[EE,abs(E2)];
            x=x-gap;
            
            
        end
        
    end
         
end

for iii=1:length(XX)
    figure(1);hold on; xline(XX(iii));
end

figure(1);hold on; plot(X,EE,'r');






















% z=(0:0.001:2);
% for jj=1:length(z)
%     E=10*exp(-10*z(jj));
%     EE(jj)=E;
%     ZZ(jj)=z(jj);
% end
% plot(ZZ,EE,'r');
