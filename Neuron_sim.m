clear variables; close all; clc

%% Model from Harish & Hansel: "Asynchronous Rate Chaos in Spiking Neuronal Circuits"

rng(1); %fix RNG seed

N = 320; %number of oscillators
K = 80; %avg # connections per node

C = zeros(N);
C(randperm(numel(C), N*K)) = 1;

tspan = [0 10];

opts = odeset('RelTol',1e-10);

% C = randi(2,N) - 1; %connections
J0 = 10; %Connection strength
I0 = 2; %Constant input
J = -(J0/sqrt(K))*C;
I = sqrt(K)*I0;
tau = 0.01;
phi = @(h) (1/2) * (1 + erf(h/sqrt(2)));
dhdt = @(h, I, J, tau) (1/tau) * (-h + I + J*phi(h));

h0 = rand(N,1);
[t, h] = ode45(@(t,h)dhdt(h,I,J,tau),tspan,h0,opts);

plot(t,h(:,1:3))

tStep = mean(diff(t));
nSteps = ceil(tspan(2)/tStep);
tN = 0:tStep:tspan(2);
tN = tN(1:end-1); %match to nSteps

h = interp1(t,h,tN);
t = tN;

save('neuron_sim_data.mat','h','t');

return;

%%
Gamma = 10;
sigma = 3;
c = 0;
log_thresh = 0.5;
log_k = 10;
log_amp = 1;
f = @(x,log_thresh,log_k,log_amp) sum(log_amp./(1 + exp(-log_k*(x-log_thresh))));

% Visualize f
% x = linspace(-2,3,100);
% y = f(x,log_thresh,log_k,log_amp);
% plot(x,y)

dxdt = @(x,Gamma,sigma,log_thresh,log_k,log_amp) -Gamma*x + sigma*f(x,log_thresh,log_k,log_amp) - c;
x0 = rand(N,1);

[t, x] = ode45(@(t,x)dxdt(x,Gamma,sigma,log_thresh,log_k,log_amp),tspan,x0,opts);

plot(t,x)

%% Plot Poincare Section
cut1 = -10.5;

figure
for j = 1:nSteps-1
    if h(j,1) <= cut1 && h(j+1,1) >= cut1
        plot(h(j,2),h(j,3),'r.');
        hold on
    elseif h(j,1) >= cut1 && h(j+1,1) <= cut1
        plot(h(j,2),h(j,3),'b.');
        hold on
    end
end