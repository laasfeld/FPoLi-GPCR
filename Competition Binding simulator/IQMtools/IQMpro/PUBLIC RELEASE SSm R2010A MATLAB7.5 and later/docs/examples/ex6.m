function J=ex6(u,tfinal,xx,rho)

%Set integration options
tol=1e-7;
options = odeset('RelTol',tol,'AbsTol',tol*ones(1,2));
%Integration
[tout,yout] = ode15s(@dynamic_model,[0 tfinal],[0.95 0.05],[],xx,u',tfinal,rho);

J=-yout(end,2);       %Change sign since we are maximizing

%***********************************************************
%Function of the dynamic system
function dy=dynamic_model(t,y,xx,u,tfinal,rho)

dy=zeros(2,1);  %Initialize the state variables

%Linear interpolation of the control variable
T = interp1q(xx,u,t);

k1=5.35e10*exp(-9000/T);
k2=4.61e17*exp(-15000/T);

dy(1)=-k1*y(1);
dy(2)=k1*y(1)-k2*y(2);
%***********************************************************