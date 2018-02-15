clear variable; close all; clc

rng(1); %fix RNG seed

N = 128; %number of oscillators
K = 0.5; %coupling strength parameter
% omega = -1 + 2*((1:N) - 1)/(N-1); %natural frequencies of oscillators
% omega = ones(1,N);
sigma = 0.2;
omega = randn(1,N)*sigma;
psi0 = pi*rand(1,N);
% psi0_2 = psi0 + 0.1*rand(1,N);

tspan = [0 2000];

K_range = 0.05:0.05:0.7;
r_inf = zeros(size(K_range));

for j = 1:length(K_range)
    K = K_range(j);
    c = K/N;
    W = omega.';
    dPdt = @(P,W,c) W+c*sum(sin(meshgrid(P)-meshgrid(P)'),2);
    [~, P] = ode45(@(t,p)dPdt(p,W,c),tspan,psi0);
    R = (1/N) * sum(exp(sqrt(-1)*P),2);
    r_inf(j) = mean(abs(R(end-50:end)));
end

figure
plot(K_range,r_inf)

%%
tspan = [0 2000];

Kc = 0.35;
K = Kc;
c = K/N;
W = omega.';
dPdt = @(P,W,c) W+c*sum(sin(meshgrid(P)-meshgrid(P)'),2);
[T, P] = ode45(@(t,p)dPdt(p,W,c),tspan,psi0);
R = (1/N) * sum(exp(sqrt(-1)*P),2);
phi = mod(diff(P,1,2),2*pi);
figure
plot(phi(:,1),phi(:,2),'.')
figure
plot(T,abs(R))