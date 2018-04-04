clear variables; close all; clc

rho = 28;
sigma = 10;
beta = 8/3;

tspan = [0 40];

opts = odeset('RelTol',1e-10);

dxdt = @(x,sigma,rho,beta) [sigma*(x(2)-x(1)); x(1)*(rho-x(3))-x(2); x(1)*x(2)-beta*x(3)];
x0 = rand(3,1);

[t, x] = ode45(@(t,x)dxdt(x,sigma,rho,beta),tspan,x0,opts);

tStep = (1/2)*mean(diff(t));
nSteps = ceil(tspan(2)/tStep);
tN = 0:tStep:tspan(2);
tN = tN(1:end-1); %match to nSteps

x = interp1(t,x,tN,'spline');
t = tN;

save('Lorenz_sim_data.mat','x','t');