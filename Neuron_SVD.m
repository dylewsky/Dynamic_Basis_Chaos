clear variables; close all; clc

load neuron_sim_data.mat;

windows = 10.^(2:0.5:3.5);
stepSize = 50;

maxSlide = floor((size(h,1) - min(windows))/stepSize);
SVD_res = cell(length(windows),maxSlide);

r = 5;

for n = 1:length(windows)
    wSteps = windows(n);
    nSlide = floor((size(h,1) - wSteps)/stepSize);
    
    figure
    for k = 1:nSlide
        thisWind = (k-1)*stepSize + 1 :(k-1)*stepSize + wSteps;
        hWind = h(thisWind,:);
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