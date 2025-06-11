clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))

%% Set inputs
split=1; % Split 1 is mild TBI set
numClusters= 6;
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
TR = 171;
load CogDemo_174sub.mat diag % vector specifying control vs. TBI
%% Separate Partition by Group

load editTP.mat unnamed1 % vector specifying TP of TBI subjects
scanInd0 = repelem(diag,TR); % for HC
scanInd1 = repelem(unnamed1(:,1),TR);% for TBI TP1
scanInd2 = repelem(unnamed1(:,2),TR);% for TBI TP2
scanInd3 = repelem(unnamed1(:,3),TR);% for TBI TP3
scanInd4 = repelem(unnamed1(:,4),TR);% for TBI TP4
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

clusterNames = {'DMN+/FPN-', 'DMN-/VATSOM+', 'DMN+/VATSOM-', 'DMN-/FPN+', 'VIS+', 'VIS-'};
%{'DMN+', 'DMN-', 'VIS-', 'VIS+'};
% Partitions of time series
HCPartition= partition(scanInd0 == 0);
T1Partition= partition(scanInd1 == 1);
T2Partition= partition(scanInd2 == 1);
T3Partition= partition(scanInd3 == 1);
T4Partition= partition(scanInd4 == 1);

% nsubj X TR vectors
HCsubjs= 39;
HCsubjInd= repelem(1:HCsubjs, TR)';

T1subjs = 30;
T1subjInd= repelem(1:T1subjs, TR)';

T2subjs = 45;
T2subjInd= repelem(1:T2subjs, TR)';

T3subjs = 34;
T3subjInd= repelem(1:T3subjs, TR)';

T4subjs = 26;
T4subjInd= repelem(1:T4subjs, TR)';
%% transition probabilities (no persist)
[HCTransitionProbability2D,HCTransitionProbabilityMats] = GET_TRANS_PROBS(HCPartition, HCsubjInd);
[HCTransitionProbabilityNoPersist2D,HCTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(HCPartition, HCsubjInd);

[T1TransitionProbability2D,T1TransitionProbabilityMats] = GET_TRANS_PROBS(T1Partition, T1subjInd);
[T1TransitionProbabilityNoPersist2D,T1TransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(T1Partition, T1subjInd);

[T2TransitionProbability2D,T2TransitionProbabilityMats] = GET_TRANS_PROBS(T2Partition, T2subjInd);
[T2TransitionProbabilityNoPersist2D,T2TransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(T2Partition, T2subjInd);

[T3TransitionProbability2D,T3TransitionProbabilityMats] = GET_TRANS_PROBS(T3Partition, T3subjInd);
[T3TransitionProbabilityNoPersist2D,T3TransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(T3Partition, T3subjInd);

[T4TransitionProbability2D,T4TransitionProbabilityMats] = GET_TRANS_PROBS(T4Partition, T4subjInd);
[T4TransitionProbabilityNoPersist2D,T4TransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(T4Partition, T4subjInd);



grpAvgHC = squeeze(mean(HCTransitionProbabilityMats,1));
% grpAvgHC = [grpAvgHC(:,3) grpAvgHC(:,2) grpAvgHC(:,4) grpAvgHC(:,1)]
% grpAvgHC = [grpAvgHC(3,:); grpAvgHC(2,:); grpAvgHC(4,:); grpAvgHC(1,:)]
%%
grpAvgT1 = squeeze(nanmean(T1TransitionProbabilityMats,1));
% grpAvgT1 = [grpAvgT1(:,3) grpAvgT1(:,2) grpAvgT1(:,4) grpAvgT1(:,1)]
% grpAvgT1 = [grpAvgT1(3,:); grpAvgT1(2,:); grpAvgT1(4,:); grpAvgT1(1,:)]

grpAvgT2 = squeeze(nanmean(T2TransitionProbabilityMats,1)); 
% grpAvgT2 = [grpAvgT2(:,3) grpAvgT2(:,2) grpAvgT2(:,4) grpAvgT2(:,1)]
% grpAvgT2 = [grpAvgT2(3,:); grpAvgT2(2,:); grpAvgT2(4,:); grpAvgT2(1,:)]

grpAvgT3 = squeeze(nanmean(T3TransitionProbabilityMats,1));
% grpAvgT3 = [grpAvgT3(:,3) grpAvgT3(:,2) grpAvgT3(:,4) grpAvgT3(:,1)]
% grpAvgT3 = [grpAvgT3(3,:); grpAvgT3(2,:); grpAvgT3(4,:); grpAvgT3(1,:)]

grpAvgT4 = squeeze(nanmean(T4TransitionProbabilityMats,1));
% grpAvgT4 = [grpAvgT4(:,3) grpAvgT4(:,2) grpAvgT4(:,4) grpAvgT4(:,1)]
% grpAvgT4 = [grpAvgHC(3,:); grpAvgT4(2,:); grpAvgT4(4,:); grpAvgT4(1,:)]

%% permutation testing to compare transition probability matrices 
[h1,pavg1,c1,t1]=ttest2(HCTransitionProbabilityMats,T1TransitionProbabilityMats); 
fdravg1 = mafdr(reshape(pavg1,1,numClusters^2),'BHFDR',1);
fdravg1 = reshape(fdravg1,numClusters,numClusters);

[h2,pavg2,c2,t2]=ttest2(HCTransitionProbabilityMats,T2TransitionProbabilityMats); 
fdravg2 = mafdr(reshape(pavg2,1,numClusters^2),'BHFDR',1);
fdravg2 = reshape(fdravg2,numClusters,numClusters);

[h3,pavg3,c3,t3]=ttest2(HCTransitionProbabilityMats,T3TransitionProbabilityMats); 
fdravg3 = mafdr(reshape(pavg3,1,numClusters^2),'BHFDR',1);
fdravg3 = reshape(fdravg3,numClusters,numClusters);

[h4,pavg4,c4,t4]=ttest2(HCTransitionProbabilityMats,T4TransitionProbabilityMats); 
fdravg4 = mafdr(reshape(pavg4,1,numClusters^2),'BHFDR',1);
fdravg4 = reshape(fdravg4,numClusters,numClusters);

%% plot transition probabilities
maxVal = max(max([grpAvgT1, grpAvgT2])); 

% Persist TP for each group
f = figure;
subplot(1,5,1);
imagesc(grpAvgHC);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('HC');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,5,2);
imagesc(grpAvgT1);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('TBI TP1');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,5,3);
imagesc(grpAvgT2);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('TBI TP2');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,5,4);
imagesc(grpAvgT3);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('TBI TP3');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,5,5);
imagesc(grpAvgT4);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('TBI TP4');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar
clusterNames = {'DMN+/FPN-', 'DMN-/VATSOM+', 'DMN+/VATSOM-', 'DMN-/FPN+', 'VIS+', 'VIS-'};
%{'DMN+', 'DMN-', 'VIS-', 'VIS+'};
% {'DMN+/FPN-', 'DMN-/VATSOM+', 'DMN+/VATSOM-', 'DMN-/FPN+', 'VIS+', 'VIS-'};

% HC - TBI TP Matrices
f= figure;
subplot(1,4,1);
HCMinusTP1 = (grpAvgHC-grpAvgT1);
imagesc(squeeze(t1.tstat),[-2 2]); colormap('plasma'); 
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Final State');
sig_thresh = 0.05; 
[y,x] = find(fdravg1 < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(pavg1) < sig_thresh); 
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
COLOR_TICK_LABELS(true,true,numClusters);
title('HC vs. 1 week post-TBI');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,4,2);
HCMinusTP2 = (grpAvgHC-grpAvgT2);
imagesc(squeeze(t2.tstat),[-2 2]); colormap(gca, 'plasma'); 
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Final State');
sig_thresh = 0.05; 
[y,x] = find(fdravg2 < sig_thresh);
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(pavg2) < sig_thresh);
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
COLOR_TICK_LABELS(true,true,numClusters);
title('HC vs. 1 month post-TBI');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,4,3);
HCMinusTP3 = (grpAvgHC-grpAvgT3);
imagesc(squeeze(t3.tstat),[-2 2]); colormap(gca, 'plasma'); 
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Final State');
sig_thresh = 0.05; 
[y,x] = find(fdravg3 < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(pavg3) < sig_thresh); 
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
COLOR_TICK_LABELS(true,true,numClusters);
title('HC vs. 6 months post-TBI');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
% 
subplot(1,4,4);
HCMinusTP4 = (grpAvgHC-grpAvgT4);
imagesc(squeeze(t4.tstat),[-2 2]); colormap(gca, 'plasma'); 
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Final State');
sig_thresh = 0.05; 
[y,x] = find(fdravg4 < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(pavg4) < sig_thresh); 
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
h = colorbar; ylabel(h,'t-statistic: HC - TBI'); 
COLOR_TICK_LABELS(true,true,numClusters);
title('HC vs. 12 months post-TBI');
set(gca,'FontSize',10);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial'); 
f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['TPvsHCTransProbs_k',num2str(numClusters),'.pdf']),'pdf');

