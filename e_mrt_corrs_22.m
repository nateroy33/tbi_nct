clear all
close all

% basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))



%% Moderate TBI Dat

% clear all; close all;

load mod_mRT.mat
load correlation_metrics_moderate.mat
HCsubjs= 21;
T1subjs= 28;
T2subjs= 11;


%%
% Correlation b/w mean TE (overall) and MRT
[R3,P3]= corr(avg_E_z, MRT2, 'rows', 'complete')
cvec4 = [repmat([255 215 0],HCsubjs,1); repmat([1 0 0],T1subjs,1); repmat([0 0 1],T2subjs,1)];
%HC gold, T1 red, T2 blue
figure; scatter(avg_E_z, MRT2,[],cvec4); lsline
title('AVG E vs MRT (Mod)');
text(4, -3,['R = ',num2str(R3);'p = ',num2str(P3)])
xlabel('AVG E');
ylabel('MRT');
