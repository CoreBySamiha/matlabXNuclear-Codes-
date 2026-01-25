function [t,NP,ND,ND2]=rk4_2(dN1,dN2,dN3,to,tf,NPo,NDo,ND2o,m)
%4th order Runge Kutta method for sloving ordinary differential equation,

h=(tf-to)/m;
t=to:h:tf;
NP=zeros(1,m+1);
ND=zeros(1,m+1);
ND2=zeros(1,m+1);
NP(1)=NPo;
ND(1)=NDo;
ND2(1)=ND2o;
for i=1:m
    k1=dN1(t(i),NP(i),ND(i),ND2(i));
    l1=dN2(t(i),NP(i),ND(i),ND2(i));
    m1=dN3(t(i),NP(i),ND(i),ND2(i));
    k2=dN1(t(i)+0.5*h,NP(i)+0.5*k1,ND(i)+0.5*l1,ND2(i)+0.5*m1);
    l2=dN2(t(i)+0.5*h,NP(i)+0.5*k1,ND(i)+0.5*l1,ND2(i)+0.5*m1);
    m2=dN3(t(i)+0.5*h,NP(i)+0.5*k1,ND(i)+0.5*l1,ND2(i)+0.5*m1);
    k3=dN1((t(i)+0.5*h),(NP(i)+0.5*k2),(ND(i)+0.5*l2),(ND2(i)+0.5*m2));
    l3=dN2((t(i)+0.5*h),(NP(i)+0.5*k2),(ND(i)+0.5*l2),(ND2(i)+0.5*m2));
    m3=dN3((t(i)+0.5*h),(NP(i)+0.5*k2),(ND(i)+0.5*l2),(ND2(i)+0.5*m2));
    k4=dN1((t(i)+h),(NP(i)+k3),(ND(i)+l3),(ND2(i)+m3));
    l4=dN2((t(i)+h),(NP(i)+k3),(ND(i)+l3),(ND2(i)+m3));
    m4=dN3((t(i)+h),(NP(i)+k3),(ND(i)+l3),(ND2(i)+m3));
    
    NP(i+1)=NP(i)+(h/6)*(k1+2*k2+2*k3+k4);
    ND(i+1)=ND(i)+(h/6)*(l1+2*l2+2*l3+l4);
    ND2(i+1)=ND2(i)+(h/6)*(m1+2*m2+2*m3+m4);

end
end