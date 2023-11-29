clear all
clear mex
%========================= PROBLEM SPECIFICATIONS =========================
rho=10;        %Number of nodes for the control profile
tfinal=30;     %Total operation time (put right units)

xl=320;         %Lower bound for the control
xu=350;         %Upper bound for the control

problem.f='ex6';

%Set the bounds for the optimization problem
problem.x_L=xl*ones(1,rho);
problem.x_U=xu*ones(1,rho);

opts.maxeval=2000;      %Maximum number of evaluations   

opts.local.solver=0;
opts.local.finish='fminsearch';

%Set the nodes (equally spaced)
xx=transpose(linspace(0,tfinal,rho));
%========================= END OF PROBLEM SPECIFICATIONS ==================

Results=ess_kernel(problem,opts,tfinal,xx,rho);

%Plot optimal control profile
plot(linspace(0,tfinal,rho),Results.xbest)
xlabel('Process Time (min)')
ylabel('Control profile (K) ')
axis([0 tfinal xl xu])
