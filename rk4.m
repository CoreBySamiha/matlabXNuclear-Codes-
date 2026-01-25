function [t,NP,ND]=rk4(dN1,dN2,to,tf,NPo,NDo,m)
%4th order Runge Kutta method for sloving ordinary differential equation,
%it can solve system of two differential equations
%where, dN1,dN2 are slope fucntions entered as strings
%to,tf left and right end points (here time interval)
%NPo,NDo are initial conditions
%m is the number of steps
h=(tf-to)/m;
t=to:h:tf;
NP=zeros(1,m+1);
ND=zeros(1,m+1);
NP(1)=NPo;
ND(1)=NDo;
for i=1:m
    k1=dN1(t(i),NP(i),ND(i));
    l1=dN2(t(i),NP(i),ND(i));
    k2=dN1(t(i)+0.5*h,NP(i)+0.5*k1,ND(i)+0.5*l1);
    l2=dN2(t(i)+0.5*h,NP(i)+0.5*k1,ND(i)+0.5*l1);
    k3=dN1((t(i)+0.5*h),(NP(i)+0.5*k2),(ND(i)+0.5*l2));
    l3=dN2((t(i)+0.5*h),(NP(i)+0.5*k2),(ND(i)+0.5*l2));
    k4=dN1((t(i)+h),(NP(i)+k3),(ND(i)+l3));
    l4=dN2((t(i)+h),(NP(i)+k3),(ND(i)+l3));
    
    NP(i+1)=NP(i)+(h/6)*(k1+2*k2+2*k3+k4);
    ND(i+1)=ND(i)+(h/6)*(l1+2*l2+2*l3+l4);
end
end