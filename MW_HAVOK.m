clear variables; close all; clc

dataLabel = 'Neuron';
inFile = [dataLabel '_sindy_input.mat'];
load(inFile);

nDelay = 1;
delaySteps = 7;

wind_no = 3; %which window length to use


x = V_full_discr_all{wind_no};
t = t_discr_all{wind_no};

x = x(:,1); %just keep 1st mode time series (following example of HAVOK on Lorenz)

dt = t(2)-t(1);

%% EIGEN-TIME DELAY COORDINATES
stackmax = 16; % Number of shift-stacked rows
r=16; % Rank of HAVOK Model
H = zeros(stackmax,size(x,1)-stackmax);
for k=1:stackmax
    H(k,:) = x(k:end-stackmax-1+k,1);
end

[U,S,V] = svd(H,'econ'); % Eigen delay coordinates

%% COMPUTE DERIVATIVES (4TH ORDER CENTRAL DIFFERENCE)
dV = zeros(length(V)-5,r);
for i=3:length(V)-3
    for k=1:r
        dV(i-2,k) = (1/(12 * dt)) * (-V(i+2,k)+8 * V(i+1,k)-8 * V(i-1,k)+V(i-2,k));
    end
end
% trim first and last two that are lost in derivative
V = V(3:end-3,1:r);

%% BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
Xi = V\dV;
A = Xi(1:r-1,1:r-1)';
B = Xi(end,1:r-1)';

%% FIGURES
figure('Position',[100 75 600 900])
subplot(3,4,1:3)
imagesc(A)
caxis([min(min(Xi)) max(max(Xi))]);
title('A Matrix')

subplot(3,4,4)
imagesc(B)
colorbar
caxis([min(min(Xi)) max(max(Xi))]);
title('B Vector')

subplot(3,4,5:8)
plot(t,x)
title('Input Data (Mode 1 of MW-SVD)')
xlabel('t')
ylabel('x')

subplot(3,4,9:12)
plot(t(3:end-stackmax-3),V(:,end))
title('Control Input (v_r)')
xlabel('t')
ylabel('v_r')