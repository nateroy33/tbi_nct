clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 6;
split=3; % i am using split to denote different sets of data 1= mild TBI, 2= mod TBI initial, 3= mod TBI w/ additional subjs
savedir = fullfile(basedir,'results','example2');
mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
%load MILDcentroids.mat
% load MODcentroids.mat
% load Dcentroids.mat
% load A4Fcentroids.mat
% load m4.mat
% load OCmild.mat
% load DcentroidsMOD.mat
% load DcentroidsK6.mat
TR = 420;

%%
overallNames = {'SUB-', 'SUB+', 'DMN+/SOM-', 'DMN-/SOM+', 'VIS+', 'VIS-'};
%{'DMN+', 'SUB-', 'VIS+', 'VIS-'}
%{'DMN+/FPN-', 'DMN-/VATSOM+', 'DMN+/VATSOM-', 'DMN-/FPN+', 'VIS+', 'VIS-'}
%{'DMN+', 'DMN-', 'VIS-', 'VIS+'}%, 'MS-3'};
% centroids= centroids;
[nparc,numClusters] = size(centroids);
% 
[~,~,~,net9angle] = NAME_CLUSTERS_ANGLE(centroids);
% 
YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','SUB+','CB+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-','SUB-','CB-'};

%% plot

clusterColors = GET_CLUSTER_COLORS(numClusters);

YeoColors = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
YeoColors = [YeoColors;YeoColors];

f = figure;
imagesc(net9angle); ax = gca;
colormap('plasma')
set(ax,'xaxisLocation','top')
xticks(1:18); xticklabels(YeoNetNames);
xtickangle(90);
for K = 1:18
	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	YeoColors(K,:), ax.XTickLabel{K});
end
yticks(1:numClusters); COLOR_TICK_LABELS(false,true,numClusters);
set(ax,'FontSize',8);

%%

[~,~,net9angle_Up,net9angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids);
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN', 'SUB', 'CB'};
numNets = numel(YeoNetNames);




%% make radial plots

clusterColors = GET_CLUSTER_COLORS(numClusters);

clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{10} = '';

f=figure;
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net9angle_Up(K,:) net9angle_Up(K,1)],'k');
    polarplot(netAngle,[net9angle_Down(K,:) net9angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rlim([0.0 0.8]);
    rticks([0.2 0.4 0.8]); rticklabels({'','0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];
