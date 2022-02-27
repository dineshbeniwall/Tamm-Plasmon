function outm=m1(lmd,th,d)
    n2=1.5;
    n1=2.8;
    np=[n1,n2];
    a=150e-9;
    
    k=2*pi/lmd;
    c=3e8;
%     Gold
    lmdp=168e-9;
    lmdc=8934e-9;
%     Silver
%     lmdp=145.41e-9;
%     lmdc=17614e-9;
    wp=2*pi*c/lmdp;
    wc=2*pi*c/lmdc;
    w=k*c; 
    nm=sqrt(1-(wp^2/(w^2-1i*w*wc)));
    ky1=k*np(1)*sind(th);
    ky2=k*nm*sind(th);
    k1z=sqrt((np(1)*w/c)^2-ky1^2);
    k2z=sqrt((nm*w/c)^2-ky2^2);
    A1= exp(1i*k1z*a)*(cos(k2z*d)+1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*d));
    B1= exp(-1i*k1z*a)*(1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*d));
    D1= exp(-1i*k1z*a)*(cos(k2z*d)-1i/2*(k2z/k1z+k1z/k2z)*sin(k2z*d));
    C1= exp(1i*k1z*a)*(-1i/2*(k2z/k1z-k1z/k2z)*sin(k2z*d));
    outm=[A1,B1;C1,D1];
end
