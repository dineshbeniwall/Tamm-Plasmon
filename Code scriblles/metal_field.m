clc;clear all;close all;

n2=1.5;
n1=2.8;
np=[n1,n2];
a=180e-9;
b=80e-9;
d=50e-9;
lmd=600e-9;
k=2*pi/lmd;
c=3e8;

w=k*c; 
v=a+b;
% at angle theta ky=k*sin(th)
th=0;

ky1=k*np(1)*sind(th);  
ky2=k*np(2)*sind(th);
k1z=sqrt((np(1)*w/c)^2-ky1^2);
k2z=sqrt((np(2)*w/c)^2-ky2^2);
A= exp(1i*k1z*a)*(cos(k2z*b)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
B= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
C= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*b));
D= exp(-1i*k1z*a)*(cos(k2z*b)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*b));
  
K=1/v*acos((A+D)/2);
a0=B;
b0=exp(1i*K*v)-A;
outm=m1(lmd,th,d);
% numbere of layers n
n=10;

z=(0:0.01e-9:n*v);
zz=0;
q=0;
gap=0.1e-9;


XX=[];
Z=[];
EE=[];
outn=[a0,b0];
    
%     
% for ii=1:n
%     outn=[D,-B;-C,A]^q*[a0;b0];
%     cndn=[exp(1i*k1z*a),exp(-1i*k1z*a);k1z/k2z*exp(1i*k1z*a),-k1z/k2z*exp(-1i*k1z*a)]*[-exp(-1i*k2z*a),-exp(-1i*k2z*a);-exp(1i*k2z*a),exp(1i*k2z*a)]*outn;
%     
%     for jj=1:v/gap
%         if jj<=b/gap
%             E2=abs(cndn(1)*exp(-1i*k2z*(zz-q*v))+cndn(2)*exp(1i*k2z*(zz-q*v)));
%             Z=[Z,zz];
%             EE=[EE,E2];
%             zz=zz+gap;
%         else          
%            E1=abs(outn(1)*exp(-1i*k1z*(zz-q*v))+outn(2)*exp(1i*k1z*(zz-q*v)));
%             Z=[Z,zz];
%             EE=[EE,E1];
%             zz=zz+gap;
%         end
%         if jj==b/gap
%             XX=[XX,zz];
%         elseif jj==v/gap
%             XX=[XX,zz];
%         end
%     end
%     q=q+1;
%     
%     
% end
% for iii=1:length(XX)
%     figure(1);hold on; xline(XX(iii));
% end
% 
% figure(1);hold on; plot(Z,EE,'r');


for ii=1:n
    
    for jj=1:v/gap
        if jj<=b/gap
            E2=abs(((outn(1)*exp(-1i*k2z*(zz-q*v-a))+outn(2)*exp(1i*k2z*(zz-q*v-a)))*exp(1i*K*(zz-q*v-a)))*exp(1i*K*zz));
            Z=[Z,zz];
            EE=[EE,E2];
            zz=zz+gap;
        else          
            E1=abs(((outn(1)*exp(-1i*k1z*(zz-q*v))+outn(2)*exp(1i*k1z*(zz-q*v)))*exp(1i*K*(zz-q*v)))*exp(1i*K*zz));
            Z=[Z,zz];
            EE=[EE,E1];
            zz=zz+gap;
        end
        if jj==b/gap
            XX=[XX,zz];
        elseif jj==v/gap
            XX=[XX,zz];
        end
    end
    q=q+1;
    
    
end
for iii=1:length(XX)
    figure(1);hold on; xline(XX(iii));
end

figure(1);hold on; plot(Z,EE,'r');


%%%long way method%%%%
% for jj=1:length(z)
%     
%     outn=[D,-B;-C,A]^q*[a0;b0];
%     cndn=[exp(1i*k1z*a),exp(-1i*k1z*a);k1z/k2z*exp(1i*k1z*a),-k1z/k2z*exp(-1i*k1z*a)]*[-exp(-1i*k2z*a),-exp(-1i*k2z*a);-exp(1i*k2z*a),exp(1i*k2z*a)]*outn;
%     E1=abs(outn(1)*exp(-1i*k1z*(z(jj)-q*v))+outn(2)*exp(1i*k1z*(z(jj)-q*v)));
%     E2=abs(cndn(1)*exp(-1i*k2z*(z(jj)-q*v))+cndn(2)*exp(1i*k2z*(z(jj)-q*v)));
%     
%     Z(jj)=z(jj);
%     EE(jj)=E; 
%     if mod(jj,26) == 0
%         q=q+1;
%         q
%     end
%     
% end
% figure(1);hold on; plot(Z,EE,'r');





