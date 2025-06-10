
clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
E= load ('TS86.mat'); % make sure this file matches the split you want to run below
savedir = fullfile(basedir,'results','elbowr2');mkdir(savedir);		% set save directory
%% inputs

concTS= E.TS86; %make sure this matches the split
N = size(concTS,1); % number of observations

split = 1; % i am using split to denote different data sets
k_rng = 2:13;
VarianceExplained = zeros(length(k_rng),1);

for numClusters = k_rng
	disp(['K = ',num2str(numClusters)])
	load(fullfile(basedir,['/repkmeans/kmeans',num2str(split),'k_',num2str(numClusters),'.mat']));
	kClusterCentroids = GET_CENTROIDS(concTS,partition,numClusters);
	VarianceExplained(numClusters - 1) = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

save(fullfile(savedir,['VarianceExplained',num2str(split),'.mat']),'VarianceExplained','k_rng');

% Fig S2a-b
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['ElbowPlot',num2str(split),'.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['GainInVarianceExplained',num2str(split),'.pdf']));