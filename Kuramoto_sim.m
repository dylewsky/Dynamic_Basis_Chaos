clear variables; close all; clc

rng(1); %fix RNG seed

N = 16; %number of oscillators
% K = 0.2; %coupling strength parameter
% omega = -1 + 2*((1:N) - 1)/(N-1); %natural frequencies of oscillators
% omega = ones(1,N);
sigma = 0.2;
% omega = randn(1,N)*sigma;
omega = linspace(-0.15,0.15,N);
% psi0 = pi*rand(1,N);
psi0 = zeros(1,N);
% psi0_2 = psi0 + 0.1*rand(1,N);

opts = odeset('RelTol',1e-7);


tspan = [0 2000];

K_range = [0.05:0.05:0.15, 0.16:0.01:0.25, 0.30:0.05:0.70];
r_inf = zeros(size(K_range));

for j = 1:length(K_range)
    K = K_range(j);
    c = K/N;
%     W = omega.';
%     dPdt = @(P,W,c) W+c*sum(sin(meshgrid(P)-meshgrid(P)'),2);
%     [~, P] = ode45(@(t,p)dPdt(p,W,c),tspan,psi0);
    [~, P] = ode45(@(t,p) Kuramoto_RHS(p,omega,K),tspan,psi0,opts);
    R = (1/N) * sum(exp(sqrt(-1)*P),2);
    r_inf(j) = mean(abs(R(end-50:end)));
end

figure
plot(K_range,r_inf)
xlabel('K');
ylabel('|R|');

%%
tspan = [0 10000];
% tspan = 0:0.0001:1000;
% Kc = 0.35;
K = 0.175;
% K = 0.05;
hold on
plot([K K],[0 1],'k:')
hold off
c = K/N;
% W = omega.';
% dPdt = @(P,W,c) W+c*sum(sin(meshgrid(P)-meshgrid(P)'),2);
% [T, P] = ode45(@(t,p)dPdt(p,W,c),tspan,psi0);
[T, P] = ode45(@(t,p) Kuramoto_RHS(p,omega,K),tspan,psi0,opts);
R = (1/N) * sum(exp(sqrt(-1)*P),2);
phi = diff(P,1,2);
% phi1_slice = abs(diff(phi(:,1))) > (2*pi - 2);
np = 1;
ps = zeros(ceil(max(phi(:,1))/(2*pi)),2);
for j = 1:length(T)
    if phi(j,1) >= 2*pi*np
        ps(np,1) = phi(j,2);
        ps(np,2) = phi(j,3);
        np = np+1;
    end
end
ps = mod(ps,2*pi);
figure
subplot(1,2,1)
plot(T,P)
axis square
subplot(1,2,2)
% plot(mod(phi(:,2),2*pi),mod(phi(:,3),2*pi),'k.')
axis square
hold on
% plot(phi(phi1_slice,2),phi(phi1_slice,3),'r.')
plot(ps(:,1),ps(:,2),'r.')
% xlim([0 2*pi]);
% ylim([0 2*pi]);
% figure
% plot(T,abs(R))

% figure
% plot(T,cos(phi))
% xlim([1000 1100])

tStep = mean(diff(T))/2;
nSteps = ceil(tspan(2)/tStep);
tN = 0:tStep:tspan(2);
tN = tN(1:end); %match to nSteps

h = interp1(T,cos(phi),tN); %enforce evenly spaced time steps
t = tN;

% figure
% plot(t,h)
% xlim([1000 1100])

save('Kuramoto_sim_data.mat','h','t');