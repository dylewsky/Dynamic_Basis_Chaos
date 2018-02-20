clear variables; close all; clc

load neuron_sim_data.mat;
h = h.';
h_const = mean(h,2);
h = h-repmat(h_const,1,length(t));

windows = floor(10.^(2:0.25:3));
stepSize = 50;

maxSlide = floor((size(h,2) - min(windows))/stepSize);
SVD_res = cell(length(windows),maxSlide);

r = 5;

for n = 1:length(windows)
    wSteps = windows(n);
    nSlide = floor((size(h,2) - wSteps)/stepSize);
    disp(['Running n = ' num2str(n)])
    figure
    for k = 1:nSlide
        thisWind = (k-1)*stepSize + 1 :(k-1)*stepSize + wSteps;
        hWind = h(:,thisWind);
        tWind = t(thisWind);
        [Uw,Sw,Vw] = svd(hWind,'econ');
        SVD_res{n,k}.U = Uw(:,1:r);
        SVD_res{n,k}.V = Vw(:,1:r);
        SVD_res{n,k}.S = diag(Sw);
        semilogy(diag(Sw));
        hold on
    end
    title(['Spectra for ' num2str(wSteps*(t(2)-t(1))) ' Second Window'])
end
save('SVD_res.mat','SVD_res','windows','stepSize','r');

%% Compare SVD Spectra
addpath(genpath('kakearney-boundedline'));

allMeans = zeros(length(windows),min(windows)); %truncate all spectra to length of shortest
allStds = zeros(size(allMeans));
b = zeros(min(windows),1,length(windows));

for n = 1:length(windows)
    wSteps = windows(n);
    nSlide = floor((size(h,1) - wSteps)/stepSize);
    disp(nSlide)
    all_S = zeros(wSteps,nSlide);

    for k = 1:nSlide
        all_S(:,k) = SVD_res{n,k}.S/sum(SVD_res{n,k}.S); %normalize
    end
    
    mean_S = mean(all_S,2);
    std_S = std(all_S,0,2);
    allMeans(n,:) = mean_S(1:min(windows));
    allStds(n,:) = std_S(1:min(windows));
    b(:,1,n) = std_S(1:min(windows));
end
% allStds = reshape(allStds,length(windows),1,min(windows)); %add singleton dimension
figure('Position',[200 200 1000 400])
subplot(1,2,1)
[h1, hp] = boundedline(1:min(windows),allMeans,b,'o');
legend(string(windows));
title('SVD Spectra by Window Size')
xlim([1 20]);
subplot(1,2,2)
[h1L, hpL] = boundedline(1:min(windows),allMeans,b,'o');
set(gca,'YScale','log')
legend(string(windows));
title('SVD Spectra by Window Size (Log Scale)')
xlim([1 20]);

%% Moving Window SVD Reconstruction

% NEED TO IMPLEMENT AN ITERATIVE MATCHING SCHEME FOR WINDOW MODES
% Check if Mode 1 from Step 1 to Step 2 is epsilon-consistent in its
% direction, or if it needs to be flipped (and the corresponding V column
% flipped too)

% ALSO: currently h_recon has a bunch of NaNs on the end; this needs to be
% addressed

h_recons = cell(length(windows),1); 
V_full_all = cell(length(windows),1);
allModes = cell(length(windows),1);
windMids_all = cell(length(windows),1);

for n = 1:length(windows)
    wSteps = windows(n);
    nSlide = floor((size(h,1) - wSteps)/stepSize);
    
    h_recon = zeros(size(h));
    V_full = zeros(r,length(t));
    wModes = zeros(nSlide,r,size(h,1)); %window step #, mode #, mode vector
    wSVs = zeros(nSlide,r); %singular values over time
    
    wCount = zeros(size(t)); %count # windows contributing to each step
    windMids = zeros(nSlide,1);
    disp(['Running n = ' num2str(n)])
    
    for k = 1:nSlide
        thisWind = (k-1)*stepSize + 1 :(k-1)*stepSize + wSteps;
        windMid = (k-1)*stepSize + floor(wSteps/2);
        windMids(k) = windMid;
        V_wind = SVD_res{n,k}.V(:,1:r);
        U_wind = SVD_res{n,k}.U(:,1:r);
        S_wind = SVD_res{n,k}.S(1:r);
        V_full(:,thisWind) = V_full(:,thisWind) + V_wind.';
        h_recon(:,thisWind) = h_recon(:,thisWind) + U_wind * diag(S_wind) * V_wind.';
        wCount(thisWind) = wCount(thisWind) + 1;
        wModes(k,:,:) = U_wind.';
        wSVs(k,:) = S_wind;
    end
    h_recon = h_recon./repmat(wCount,size(h,1),1);
    V_full = V_full./repmat(wCount,r,1);
    
    h_recons{n} = h_recon;
    V_full_all{n} = V_full;
    allModes{n} = wModes;
    windMids_all{n} = windMids;
    
    figure
    subplot(2,2,1)
    plot(t,h_recon(1:r,:))
    title(['Reconstruction: ' num2str(wSteps)])
    subplot(2,2,2)
    plot(t,V_full)
    title('Mode Proj. Time Series')
    subplot(2,2,3)
    % plot top r modes' first coordinates over time
    plot(windMids*(t(2)-t(1)), squeeze(wModes(:,1:r,1)))
    title('Modes over Time')
    subplot(2,2,4)
    plot(windMids*(t(2)-t(1)),wSVs)
    title('Singular Values over Time');
end
