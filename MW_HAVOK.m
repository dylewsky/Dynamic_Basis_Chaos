clear variables; close all; clc

dataLabel = 'Neuron';
inFile = [dataLabel '_sindy_input.mat'];
load(inFile);

nDelay = 1;
delaySteps = 7;

wind_no = 3; %which window length to use

x = V_full_discr_all{wind_no};
t = t_discr_all{wind_no};
