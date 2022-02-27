close all; clear all; clc;




x=0;
E=[];
X=[];

count=4;
[Eout,Earray,Xarray,xfinal] = field(count,E,X,x);
figure(1);hold on;plot(Xarray,Earray,'r');

function [Eout,Earray,Xarray,xfinal]=field(count,E,X,x)
    n1=2.8;
    a=200e-9;
    n2=1.5;
    b=80e-9;
    n0=1;

    lmd=700e-9;
    c=3e8;
    w=2*pi*c/lmd;

    ep0=n0^2;
    ep1=n1^2;
    ep2=n2^2;
    k0=sqrt(w^2*ep0);
    k1=sqrt(w^2*ep1);
    k2=sqrt(w^2*ep2);

    th0=60;
    th1=asind(n0*sind(th0)/n1);
    th2=asind(n1*sind(th1)/n2);
    step=1e-9;
    na=a/step;
    nb=b/step;
    
    if count==0
        
        E0=ep0*sind(th0)*cosd(th1)+ep1*sind(th1)*cosd(th0);
        E1=2*ep0*sind(th0)*cosd(th0);
        E00=ep0*sind(th0)*cosd(th1)-ep1*sind(th1)*cosd(th0);
        EX0=E0*exp(-1i*k0*cosd(th0)*x)+E00*exp(1i*k0*cosd(th0)*x);
        
        X=[X,x];
        E=[E,abs(EX0)];
        Eout=E1;
        Earray=E;
        Xarray=X;
        xfinal=x;
        return;
    elseif mod(count,2)==0
        [Eout,Earray,Xarray,xfinal]=field(count-1,E,X,x);
        E=Earray;
        X=Xarray;
        x=xfinal;
        E2=Eout;
        E1=E2*(1+ep2)/(cosd(th1)/cosd(th2)+(ep1*sind(th1))/(ep2*sind(th2)));
        E22=E2*(1-(ep1*sind(th1)*(1+ep2))/(ep2*sind(th2)*(cosd(th1)/cosd(th2)+(ep1*sind(th1)/(ep2*sind(th2))))));
        for ii=1:nb
            x=x+step;
            count
            EX2=E2*exp(-1i*k2*cosd(th2)*x-(count*(a+b))-a)+E22*exp(1i*k2*cosd(th2)*x-(count*(a+b))-a);
            X=[X,x];
            E=[E,abs(EX2)];
        end
        
        Eout=E1;
        Earray=E;
        Xarray=X;
        xfinal=x;
        return;
    else
        [Eout,Earray,Xarray,xfinal]=field(count-1,E,X,x);
        E=Earray;
        X=Xarray;
        x=xfinal;
        E1=Eout;
        E2=E1*(1+ep1)/(cosd(th2)/cosd(th1)+(ep2*sind(th2))/(ep1*sind(th1)));
        E11=E1*(1-(ep2*sind(th2)*(1+ep1))/(ep1*sind(th1)*(cosd(th2)/cosd(th1)+(ep2*sind(th2)/(ep1*sind(th1))))));
        for ii=1:na
            x=x+step;
            EX1=E1*exp(-1i*k1*cosd(th1)*x-(count*(a+b)))+E11*exp(1i*k1*cosd(th1)*x-(count*(a+b)));
            X=[X,x];
            E=[E,abs(EX1)];
            
        end
        
        Eout=E2;
        Earray=E;
        Xarray=X;
        xfinal=x;
        return;
    end
    
end


% %layer 0
% E0=ep0*sind(th0)*cosd(th1)+ep1*sind(th1)*cosd(th0);
% E1=2*ep0*sind(th0)*cosd(th0);
% E00=ep0*sind(th0)*cosd(th1)-ep1*sind(th1)*cosd(th0);
% 
% %layer 1
% E1=E1;
% E2=E1*(1+ep1)/(cosd(th2)/cosd(th1)+(ep2*sind(th2))/(ep1*sind(th1)));
% E11=E1*(1-(ep2*sind(th2)*(1+ep1))/(ep1*sind(th1)*(cosd(th2)/cosd(th1)+(ep2*sind(th2)/(ep1*sind(th1))))));
% 
% %layer2
% E2=E2;
% E1=E2*(1+ep2)/(cosd(th1)/cosd(th2)+(ep1*sind(th1))/(ep2*sind(th2)));
% E22=E2*(1-(ep1*sind(th1)*(1+ep2))/(ep2*sind(th2)*(cosd(th1)/cosd(th2)+(ep1*sind(th1)/(ep2*sind(th2))))));
% 
% %layer 3 will be same as layer 1 
% %layer 4 will be same as layer 2
% 
% EX0=E0*exp(-1i*k0*cosd(th0)*x)+E00*exp(1i*k0*cosd(th0)*x);
% EX1=E1*exp(-1i*k1*cosd(th1)*x)+E11*exp(1i*k1*cosd(th1)*x);
% EX2=E2*exp(-1i*k2*cosd(th2)*x)+E22*exp(1i*k2*cosd(th2)*x);







