%% load stuff
clear all
close all

basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))

split = 1;
numClusters =4;
TR = 171; 
nsubjs= 174;
nparc= 86;

load TS86.mat
load FODTAR.mat 

rTS86= reshape(TS86,TR,nsubjs,nparc);
rTS86= permute(rTS86,[2,1,3]);
concTS = rTS86;

HCsubjs= 39;

subjInd=repelem(1:nsubjs,TR)';

clusterNames={};
clusterNamesUp={};
clusterNamesDown={};

%% compute centroids

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'partition','centroids'); 
partition= reshape(partition,[TR,nsubjs])';
overallCentroids=centroids;

centroids=NaN(nsubjs,nparc,numClusters);

for i=1:nsubjs
    
    disp(['Generating centroids for scan ',num2str(i)]);
    centroids(i,:,:) = GET_CENTROIDS(squeeze(concTS(i,:,:)),partition(i,:),numClusters);
end
% >>> Insert shift here for positive scale validation
min_val = min(centroids(:));
epsilon = 0.001;
shift_constant = abs(min_val) + epsilon;
centroids = centroids + shift_constant;
% %% Reorder Centroids
% c1= centroids(:,:,1); %reorder centroids based on top-down vs. bottom-up
% c2= centroids(:,:,2);
% c3= centroids(:,:,3);
% c4= centroids(:,:,4);
% 
% centroids= cat(3,c2,c4,c3,c1)% 
% 
% o1= overallCentroids(:,1); %reorder centroids based on top-down vs. bottom-up
% o2= overallCentroids(:,2);
% o3= overallCentroids(:,3);
% o4= overallCentroids(:,4);
% 
% overallCentroids= cat(2,o2,o4,o3,o1)% 

%% load sc etc for E_calcs

load SCmats174.mat
load CogDemo_174sub.mat diag

HCSCmats= SC(:,:,diag==0);
HCSCavg = nanmean(HCSCmats,3);

c = 0; 

Anorm = NORMALIZE(HCSCavg,c);


Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];
x0 = overallCentroids(:,Xo_ind);
xf = overallCentroids(:,Xf_ind);
% Set T
T=1;

%% Flavor 1: Indiv SC/Indiv Centroids

E_full=NaN(nsubjs,numClusters^2); % !!!!!

for i=1:nsubjs
    sci = squeeze(SC(:,:,i));
    sci(find(eye(nparc,nparc)))= 0;
    if ~isnan(sci(2,1))
        Anorm = NORMALIZE(sci,c);
        
        %define initial and final states for each subject/scan
        x0 = squeeze(centroids(i,:,Xo_ind));
        xf = squeeze(centroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
        WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
        E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    end
    
end
load editTP.mat unnamed1 % matrix specifying TP of TBI subjects
scanInd0 = diag; % for HC
scanInd1 = unnamed1(:,1);
scanInd2 = unnamed1(:,2);
scanInd3 = unnamed1(:,3);
scanInd4 = unnamed1(:,4);

HCef= E_full(scanInd0==0,:);
HCef= [HCef(1:21,:);HCef(23:end,:)]
[h1,p1,ci1,tstat1] = ttest2(HCef,E_full(scanInd1==1,:));
[h2,p2,ci2,tstat2] = ttest2(HCef,E_full(scanInd2==1,:));
[h3,p3,ci3,tstat3] = ttest2(HCef,E_full(scanInd3==1,:));
[h4,p4,ci4,tstat4] = ttest2(HCef,E_full(scanInd4==1,:));

figure;
imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat])
allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat];
pall = [p1; p2; p3; p4];
% p1= reshape(p1,[numClusters numClusters])';
p1u= mafdr(p1,'BHFDR',1)
p1= reshape(p1,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
p2u= mafdr(p2,'BHFDR',1)
p2= reshape(p2,[numClusters numClusters])';
% p3= reshape(p3,[numClusters numClusters])';
p3u= mafdr(p3,'BHFDR',1)
p3= reshape(p3,[numClusters numClusters])';
% p4= reshape(p4,[numClusters numClusters])';
p4u= mafdr(p4,'BHFDR',1)
p4= reshape(p4,[numClusters numClusters])';

maxVal = max(allt(:)); % sync color scales
minVal = min(allt(:));
sig_thresh = 0.05;
clusterNames = {'DMN+','DMN-','VIS-','VIS+'};

subplot(1,4,1);
imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 1 week post-TBI');
set(gca,'FontSize',14);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p1 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,2);
imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 1 month post-TBI');
set(gca,'FontSize',14);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
[y,x] = find(p2 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,3);
imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 6 months post-TBI');
set(gca,'FontSize',14);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p3 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,4);
imagesc(reshape(tstat4.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 12 months post-TBI');
set(gca,'FontSize',14);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p4 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end
%% Longitudinal TBI Comparisons
[h1,p1,ci1,tstat1] = ttest2(E_full(scanInd1==1,:),E_full(scanInd2==1,:));
[h2,p2,ci2,tstat2] = ttest2(E_full(scanInd2==1,:),E_full(scanInd3==1,:));
[h3,p3,ci3,tstat3] = ttest2(E_full(scanInd3==0,:),E_full(scanInd4==1,:));


figure;
imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat])
allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat];
pall = [p1; p2; p3];
p1= reshape(p1,[numClusters numClusters])';
p2= reshape(p2,[numClusters numClusters])';
p3= reshape(p3,[numClusters numClusters])';

maxVal = max(allt(:)); % sync color scales
minVal = min(allt(:));
sig_thresh = 0.05;
clusterNames = {'DMN+','DMN-','VIS-','VIS+'};
subplot(1,4,1);
imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('1 week vs 1 month post-TBI');
set(gca,'FontSize',11);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p1 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,2);
imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('1 month vs 6 months post-TBI');
set(gca,'FontSize',11);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
[y,x] = find(p2 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,3);
imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('6 months vs 12 months post-TBI');
set(gca,'FontSize',11);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p3 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

figure;

[h1,p1,ci1,tstat1] = ttest2(E_full(scanInd1==1,:),E_full(scanInd3==1,:));
[h2,p2,ci2,tstat2] = ttest2(E_full(scanInd1==1,:),E_full(scanInd4==1,:));
[h3,p3,ci3,tstat3] = ttest2(E_full(scanInd2==1,:),E_full(scanInd4==1,:));

imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat])
allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat];
pall = [p1; p2; p3];
p1= reshape(p1,[numClusters numClusters])';
p2= reshape(p2,[numClusters numClusters])';
p3= reshape(p3,[numClusters numClusters])';

maxVal = max(allt(:)); 
minVal = min(allt(:));
sig_thresh = 0.05;
clusterNames = {'DMN+','DMN-','VIS-','VIS+'};
subplot(1,4,1);
imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('TBI TP1 vs TBI TP3');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
hold on; 
[y,x] = find(p1 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end

subplot(1,4,2);
imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('TBI TP1 vs TBI TP4');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
[y,x] = find(p2 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end
subplot(1,4,3);
imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('TBI TP2 vs TBI TP4');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([-2 2]); colorbar
[y,x] = find(p3 < sig_thresh); 
if ~isempty(y)
    text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
end
%% Flavor 2: HC Avg SC/Indiv Centroids
% E_full=NaN(nsubjs,numClusters^2);
% 
% for i=1:nsubjs
%     sci = HCSCavg;
%     sci(find(eye(nparc,nparc)))= 0;
%     if ~isnan(sci(2,1))
%         Anorm = NORMALIZE(sci,c);
%         
%         %define initial and final states for each subject/scan
%         x0 = squeeze(centroids(i,:,Xo_ind));
%         xf = squeeze(centroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
%         WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
%         E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
%     end
%     
% end
% load editTP.mat unnamed1 % vector specifying TP of TBI subjects
% scanInd0 = diag; % for HC
% scanInd1 = unnamed1(:,1);
% scanInd2 = unnamed1(:,2);
% scanInd3 = unnamed1(:,3);
% scanInd4 = unnamed1(:,4);
% 
% [h1,p1,ci1,tstat1] = ttest2(E_full(scanInd0==0,:),E_full(scanInd1==1,:));
% [h2,p2,ci2,tstat2] = ttest2(E_full(scanInd0==0,:),E_full(scanInd2==1,:));
% [h3,p3,ci3,tstat3] = ttest2(E_full(scanInd0==0,:),E_full(scanInd3==1,:));
% [h4,p4,ci4,tstat4] = ttest2(E_full(scanInd0==0,:),E_full(scanInd4==1,:));
% 
% figure;
% imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat])
% allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat];
% pall = [p1; p2; p3; p4];
% p1= reshape(p1,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
% p3= reshape(p3,[numClusters numClusters])';
% p4= reshape(p4,[numClusters numClusters])';
% 
% maxVal = max(allt(:)); % sync color scales
% minVal = min(allt(:));
% sig_thresh = 0.05;
% clusterNames = {'VIS+','VIS-','DMN+','DMN-'};
% 
% subplot(1,4,1);
% imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 1');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p1 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,2);
% imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP2');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% [y,x] = find(p2 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,3);
% imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 3');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p3 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,4);
% imagesc(reshape(tstat4.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 4');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p4 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% %% Flavor 3: Indiv SC/Grp Avg Centroids
% E_full=NaN(nsubjs,numClusters^2);
% 
% for i=1:nsubjs
%     sci = squeeze(SC(:,:,i));
%     sci(find(eye(nparc,nparc)))= 0;
%     if ~isnan(sci(2,1))
%         Anorm = NORMALIZE(sci,c);
%         
%         %define initial and final states for each subject/scan
%         x0 = overallCentroids(:,Xo_ind);
%         xf = overallCentroids(:,Xf_ind); % now each column of x0 and xf represent state transitions
%         WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
%         E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
%     end
%     
% end
% load editTP.mat unnamed1 % vector specifying TP of TBI subjects
% scanInd0 = diag; % for HC
% scanInd1 = unnamed1(:,1);
% scanInd2 = unnamed1(:,2);
% scanInd3 = unnamed1(:,3);
% scanInd4 = unnamed1(:,4);
% 
% [h1,p1,ci1,tstat1] = ttest2(E_full(scanInd0==0,:),E_full(scanInd1==1,:));
% [h2,p2,ci2,tstat2] = ttest2(E_full(scanInd0==0,:),E_full(scanInd2==1,:));
% [h3,p3,ci3,tstat3] = ttest2(E_full(scanInd0==0,:),E_full(scanInd3==1,:));
% [h4,p4,ci4,tstat4] = ttest2(E_full(scanInd0==0,:),E_full(scanInd4==1,:));
% 
% figure;
% imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat])
% allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat; tstat4.tstat];
% pall = [p1; p2; p3; p4];
% p1= reshape(p1,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
% p3= reshape(p3,[numClusters numClusters])';
% p4= reshape(p4,[numClusters numClusters])';
% 
% 
% maxVal = max(allt(:)); % sync color scales
% minVal = min(allt(:));
% sig_thresh = 0.05;
% clusterNames = {'VIS+','VIS-','DMN+','DMN-'};
% 
% subplot(1,4,1);
% imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 1');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p1 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,2);
% imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP2');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% [y,x] = find(p2 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,3);
% imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 3');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p3 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% 
% subplot(1,4,4);
% imagesc(reshape(tstat4.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('HC vs TBI TP 4');
% set(gca,'FontSize',18);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% hold on; 
% [y,x] = find(p4 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end

%% Avg E by Time Point
AEhc= nanmean(E_full(scanInd0==0,:),2);
AEhc= AEhc([1:21, 23:end]);
AEt1= nanmean(E_full(scanInd1==1,:),2);
AEt2= nanmean(E_full(scanInd2==1,:),2);
AEt3= nanmean(E_full(scanInd3==1,:),2);
AEt4= nanmean(E_full(scanInd4==1,:),2);
[h1,p1,ci1,tstat1] = ttest2(AEhc,AEt4)
%%
% [h1,p1,ci1,tstat1] = ttest2(nanmean(E_full(scanInd0==0,:),2),nanmean(E_full(scanInd1==1,:),2));
% [h2,p2,ci2,tstat2] = ttest2(nanmean(E_full(scanInd0==0,:),2),nanmean(E_full(scanInd2==1,:),2));
% [h3,p3,ci3,tstat3] = ttest2(nanmean(E_full(scanInd0==0,:),2),nanmean(E_full(scanInd3==1,:),2));
% [h4,p4,ci4,tstat4] = ttest2(nanmean(E_full(scanInd0==0,:),2),nanmean(E_full(scanInd4==1,:),2));
%%
HCsubjs= 39;
T1subjs= 30;
T2subjs= 45;
T3subjs= 34;
T4subjs= 26;

subj= [HCsubjs T1subjs T2subjs T3subjs T4subjs];
maxsubj= max(subj);

data= NaN(maxsubj,5);
data(1:38,1)= nanmean(AEhc,2);
data(1:T1subjs,2)= nanmean(AEt1,2);
data(1:T2subjs,3)= nanmean(AEt2,2);
data(1:T3subjs,4)= nanmean(AEt3,2);
data(1:T4subjs,5)= nanmean(AEt4,2);

[h,p1,ci,stats] = ttest2(data(:,1),data(:,2));
[h,p2,ci,stats] = ttest2(data(:,1),data(:,3));
[h,p3,ci,stats] = ttest2(data(:,1),data(:,4));
[h,p4,ci,stats] = ttest2(data(:,1),data(:,5));
pval= [p1 p2 p3 p4];
upval = mafdr(pval,'BHFDR',1);
comp = [pval; upval];

figure;
violinplot(data);
title('Global Transition Energy');
ylabel('Transition Energy');
xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% xlabel('post-TBI')
TickSize = 16
ax = gca
ax.FontSize = TickSize
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
%%
AEhc= nanmean(E_full(scanInd0==0,:),2)
i= 1.5*iqr(AEhc)
Q1= quantile(AEhc,.25)
Q3= quantile(AEhc,.75)
F1= Q1-i
F2= Q3+i
%%

% %% By Column
% 
% VisM= E_full(:,[2 6 10 14]); %transitions into VIS-
% VisP= E_full(:,[1 5 9 13]); %transitions into VIS+
% DmnP= E_full(:,[3 7 11 15]);%transitions into DMN+
% DmnM= E_full(:,[4 8 12 16]);%transitions into DMN-
% 
% 
% [h1,p1,ci1,tstat1] = ttest2(nanmean(VisP(scanInd0==0,:),2),nanmean(VisP(scanInd1==1,:),2));
% [h2,p2,ci2,tstat2] = ttest2(nanmean(VisP(scanInd0==0,:),2),nanmean(VisP(scanInd2==1,:),2));
% [h3,p3,ci3,tstat3] = ttest2(nanmean(VisP(scanInd0==0,:),2),nanmean(VisP(scanInd3==1,:),2));
% [h4,p4,ci4,tstat4] = ttest2(nanmean(VisP(scanInd0==0,:),2),nanmean(VisP(scanInd4==1,:),2));
%%
% data= NaN(maxsubj,5);
% data(1:HCsubjs,1)= nanmean(VisP(scanInd0==0,:),2);
% data(1:T1subjs,2)= nanmean(VisP(scanInd1==1,:),2);
% data(1:T2subjs,3)= nanmean(VisP(scanInd2==1,:),2);
% data(1:T3subjs,4)= nanmean(VisP(scanInd3==1,:),2);
% data(1:T4subjs,5)= nanmean(VisP(scanInd4==1,:),2);
% % 
% figure;
% violinplot(data);
% title('Avg E for Transitions Into VIS+ by Group');
% ylabel('Avg E');
% xticklabels({'hc'; 'tbi1'; 'tbi2'; 'tbi3'; 'tbi4'});
% 
% anova1(data)
% %% Connecting Patients
% 
%% 
load 49subjE.mat 
X= [1 2 3 4]
Y= data11 %previously organized matrix of average TE for each TBI subj at each TP
figure;
plot(X,Y,'-o')
xticks([1 2 3 4])
xticklabels([{'tbi1'} {'tbi2'} {'tbi3'} {'tbi4'}])
ylabel('AvgE')
xlabel('Timepoint')
%% Top Down vs. Bottom Up
load TdBu.mat %vector indicating top-down vs. bottom-up information for each ROI
mcentroids= mean(centroids,1)
mcentroids= squeeze(mcentroids)
mcentroids= L2MAGNITUDENORM(mcentroids)
mcentroids= abs(mcentroids)
VisPD= dot(mcentroids(:,1), TdBu)
DmnPD= dot(mcentroids(:,2), TdBu)
VisMD= dot(mcentroids(:,3), TdBu)
DmnMD= dot(mcentroids(:,4), TdBu)
% LimPD= dot(mcentroids(:,5), TdBu)
% LimMD= dot(mcentroids(:,6), TdBu)

% New State Order: DMN+, VIS-, DMN-, VIS+ (3,2,4,1): defined by ordering
% above dot products from greatest to least
% %% Change in MRT vs. Change in gTE
% Ech= data11'
% Ech= (Ech-(nanmean(Ech,'all')))/(nanstd(Ech,0,'all'))
% Ech= [vertcat((Ech(:,2)-Ech(:,1)),(Ech(:,3)-Ech(:,2)),(Ech(:,4)-Ech(:,3)))]
% 
% load mild_mRT_reorderd.mat
% load MRTch.mat
% 
% MRTch= (MRTch-(nanmean(MRTch,'all')))/(nanstd(MRTch,0,'all'))
% MRTch= [vertcat((MRTch(:,2)-MRTch(:,1)),(MRTch(:,3)-MRTch(:,2)),(MRTch(:,4)-MRTch(:,3)))]
% 
% [R,P]= corr(Ech(:,1), MRTch(:,1),'rows','complete','type','spearman')
% 
% figure;
% scatter((Ech(:,1)),(MRTch(:,1)),25,'red')
% text(1, 2,['R = ',num2str(R)]);
% text(1,1.8,['p = ',num2str(P)])
% ylabel('Change in MRT')
% xlabel('Change in AVG TE')
%% Reorder E_full
% E_hc= E_full(scanInd0==0,:)
% E_hcm= (E_hc)
% HCbu= (E_hcm(:,2)+E_hcm(:,3)+E_hcm(:,4)+E_hcm(:,5)+E_hcm(:,6)+E_hcm(:,9)+E_hcm(:,10)+E_hcm(:,11)+E_hcm(:,12)+E_hcm(:,16)+E_hcm(:,17)+E_hcm(:,18)+E_hcm(:,23)+E_hcm(:,24)+E_hcm(:,30)) %(trans into more bu state)
% HCtd= (E_hcm(:,7)+E_hcm(:,13)+E_hcm(:,14)+E_hcm(:,19)+E_hcm(:,20)+E_hcm(:,21)+E_hcm(:,25)+E_hcm(:,26)+E_hcm(:,27)+E_hcm(:,28)+E_hcm(:,31)+E_hcm(:,32)+E_hcm(:,33)+E_hcm(:,34)+E_hcm(:,35)) %(trans into more td state)
% [h,p,Ci1,stats]= ttest(HCbu,HCtd)
% %%
% E_1= E_full(scanInd1==1,:)
% E_1m= (E_1)
% T1bu= (E_1m(:,2)+E_1m(:,3)+E_1m(:,4)+E_1m(:,5)+E_1m(:,6)+E_1m(:,9)+E_1m(:,10)+E_1m(:,11)+E_1m(:,12)+E_1m(:,16)+E_1m(:,17)+E_1m(:,18)+E_1m(:,23)+E_1m(:,24)+E_1m(:,30)) %(trans into more bu state)
% T1td= (E_1m(:,7)+E_1m(:,13)+E_1m(:,14)+E_1m(:,19)+E_1m(:,20)+E_1m(:,21)+E_1m(:,25)+E_1m(:,26)+E_1m(:,27)+E_1m(:,28)+E_1m(:,31)+E_1m(:,32)+E_1m(:,33)+E_1m(:,34)+E_1m(:,35)) %(trans into more td state)
% [h,p,Ci2,stats]= ttest(T1bu,T1td)
% %%
% E_2= E_full(scanInd2==1,:)
% E_2m= (E_2)
% T2bu= (E_2m(:,2)+E_2m(:,3)+E_2m(:,4)+E_2m(:,5)+E_2m(:,6)+E_2m(:,9)+E_2m(:,10)+E_2m(:,11)+E_2m(:,12)+E_2m(:,16)+E_2m(:,17)+E_2m(:,18)+E_2m(:,23)+E_2m(:,24)+E_2m(:,30)) %(trans into more bu state)
% T2td= (E_2m(:,7)+E_2m(:,13)+E_2m(:,14)+E_2m(:,19)+E_2m(:,20)+E_2m(:,21)+E_2m(:,25)+E_2m(:,26)+E_2m(:,27)+E_2m(:,28)+E_2m(:,31)+E_2m(:,32)+E_2m(:,33)+E_2m(:,34)+E_2m(:,35)) %(trans into more td state)
% [h,p,Ci3,stats]= ttest(T2bu,T2td)
% %%
% E_3= E_full(scanInd3==1,:)
% E_3m= (E_3)
% T3bu= (E_3m(:,2)+E_3m(:,3)+E_3m(:,4)+E_3m(:,5)+E_3m(:,6)+E_3m(:,9)+E_3m(:,10)+E_3m(:,11)+E_3m(:,12)+E_3m(:,16)+E_3m(:,17)+E_3m(:,18)+E_3m(:,23)+E_3m(:,24)+E_3m(:,30)) %(trans into more bu state)
% T3td= (E_3m(:,7)+E_3m(:,13)+E_3m(:,14)+E_3m(:,19)+E_3m(:,20)+E_3m(:,21)+E_3m(:,25)+E_3m(:,26)+E_3m(:,27)+E_3m(:,28)+E_3m(:,31)+E_3m(:,32)+E_3m(:,33)+E_3m(:,34)+E_3m(:,35)) %(trans into more td state)
% [h,p,Ci4,stats]= ttest(T3bu,T3td)
% %%
% E_4= E_full(scanInd4==1,:)
% E_4m= (E_4)
% T4bu= (E_4m(:,2)+E_4m(:,3)+E_4m(:,4)+E_4m(:,5)+E_4m(:,6)+E_4m(:,9)+E_4m(:,10)+E_4m(:,11)+E_4m(:,12)+E_4m(:,16)+E_4m(:,17)+E_4m(:,18)+E_4m(:,23)+E_4m(:,24)+E_4m(:,30))%(trans into more bu state)
% T4td= (E_4m(:,7)+E_4m(:,13)+E_4m(:,14)+E_4m(:,19)+E_4m(:,20)+E_4m(:,21)+E_4m(:,25)+E_4m(:,26)+E_4m(:,27)+E_4m(:,28)+E_4m(:,31)+E_4m(:,32)+E_4m(:,33)+E_4m(:,34)+E_4m(:,35)) %(trans into more td state)
% [h,p,Ci5,stats]= ttest(T4bu,T4td)

% data22= NaN(45,10);
% data22(1:39,1)= HCbu;
% data22(1:39,2)= HCtd;
% data22(1:30,3)= T1bu;
% data22(1:30,4)= T1td;
% data22(1:45,5)= T2bu;
% data22(1:45,6)= T2td;
% data22(1:34,7)= T3bu;
% data22(1:34,8)= T3td;
% data22(1:26,9)= T4bu;
% data22(1:26,10)= T4td;
% 
% figure;
% violinplot(data22);
% title('Avg TE by Group and Transition Type');
% ylabel('Avg E');
% xticklabels({'hc bu'; 'hc td'; 'tbi1 bu'; 'tbi1 td'; 'tbi2 bu'; 'tbi2 td'; 'tbi3 bu'; 'tbi3 td'; 'tbi4 bu'; 'tbi4 td';});
%%
% HCbu= (E_hcm(:,2)+E_hcm(:,3)+E_hcm(:,4)+E_hcm(:,7)+E_hcm(:,8)+E_hcm(:,12)) %(trans into more bu state)
% HCtd= (E_hcm(:,5)+E_hcm(:,9)+E_hcm(:,10)+E_hcm(:,13)+E_hcm(:,14)+E_hcm(:,15)) %(trans into more td state)
% HCr= HCbu./HCtd

% T1bu= (E_1m(:,2)+E_1m(:,3)+E_1m(:,4)+E_1m(:,7)+E_1m(:,8)+E_1m(:,12)) %(trans into more bu state)
% T1td= (E_1m(:,5)+E_1m(:,9)+E_1m(:,10)+E_1m(:,13)+E_1m(:,14)+E_1m(:,15)) %(trans into more td state)
% Tr1= T1bu./T1td

% T2bu= (E_2m(:,2)+E_2m(:,3)+E_2m(:,4)+E_2m(:,7)+E_2m(:,8)+E_2m(:,12)) %(trans into more bu state)
% T2td= (E_2m(:,5)+E_2m(:,9)+E_2m(:,10)+E_2m(:,13)+E_2m(:,14)+E_2m(:,15)) %(trans into more td state)
% Tr2= T2bu./T2td

% T3bu= (E_3m(:,2)+E_3m(:,3)+E_3m(:,4)+E_3m(:,7)+E_3m(:,8)+E_3m(:,12)) %(trans into more bu state)
% T3td= (E_3m(:,5)+E_3m(:,9)+E_3m(:,10)+E_3m(:,13)+E_3m(:,14)+E_3m(:,15)) %(trans into more td state)
% Tr3= T3bu./T3td

% T4bu= (E_4m(:,2)+E_4m(:,3)+E_4m(:,4)+E_4m(:,7)+E_4m(:,8)+E_4m(:,12)) %(trans into more bu state)
% T4td= (E_4m(:,5)+E_4m(:,9)+E_4m(:,10)+E_4m(:,13)+E_4m(:,14)+E_4m(:,15)) %(trans into more td state)
% Tr4= T4bu./T4td

% 

%% New TdBu Calc
E_hc= HCef
E_1= E_full(scanInd1==1,:)
E_2= E_full(scanInd2==1,:)
E_3= E_full(scanInd3==1,:)
E_4= E_full(scanInd4==1,:)

HCr= (((E_hc(:,5)-E_hc(:,2))+(E_hc(:,10)-E_hc(:,7))+(E_hc(:,15)-E_hc(:,12))+(E_hc(:,9)-E_hc(:,3))+(E_hc(:,14)-E_hc(:,8))+(E_hc(:,13)-E_hc(:,4)))./6)

Tr1= (((E_1(:,5)-E_1(:,2))+(E_1(:,10)-E_1(:,7))+(E_1(:,15)-E_1(:,12))+(E_1(:,9)-E_1(:,3))+(E_1(:,14)-E_1(:,8))+(E_1(:,13)-E_1(:,4)))./6)
Tr2= (((E_2(:,5)-E_2(:,2))+(E_2(:,10)-E_2(:,7))+(E_2(:,15)-E_2(:,12))+(E_2(:,9)-E_2(:,3))+(E_2(:,14)-E_2(:,8))+(E_2(:,13)-E_2(:,4)))./6)
Tr3= (((E_3(:,5)-E_3(:,2))+(E_3(:,10)-E_3(:,7))+(E_3(:,15)-E_3(:,12))+(E_3(:,9)-E_3(:,3))+(E_3(:,14)-E_3(:,8))+(E_3(:,13)-E_3(:,4)))./6)
Tr4= (((E_4(:,5)-E_4(:,2))+(E_4(:,10)-E_4(:,7))+(E_4(:,15)-E_4(:,12))+(E_4(:,9)-E_4(:,3))+(E_4(:,14)-E_4(:,8))+(E_4(:,13)-E_4(:,4)))./6)

data33= NaN(maxsubj,5);
data33(1:38,1)= HCr;
data33(1:T1subjs,2)= Tr1;
data33(1:T2subjs,3)= Tr2;
data33(1:T3subjs,4)= Tr3;
data33(1:T4subjs,5)= Tr4;

figure;
violinplot(data33);
title('Transition Energy Asymmetry');
ylabel('Bottom Up Energy - Top Down Energy');
ylim([-2.5 2])
% xlabel('post-TBI');
TickSize = 16
ax = gca
ax.FontSize = TickSize
xticklabels({'hc'; '1w post-TBI'; '1m post-TBI'; '6m post-TBI'; '12m post-TBI'});
% text(4, 1.5,['*'],'FontSize',28)
% text(5, 1.5,['*'], 'FontSize',28)
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

%%
[h,p1,Ci,stats]= ttest2(HCr,Tr1)
[h,p2,Ci,stats]= ttest2(HCr,Tr2)
[h,p3,Ci,stats]= ttest2(HCr,Tr3)
[h,p4,Ci,stats]= ttest2(HCr,Tr4)
pval= [p1 p2 p3 p4]
upval= mafdr(pval,'BHFDR',1)
%% mRT vs. TdBu
load mild_mRT_reorderd.mat
mRT= mRT'
TdBu1= vertcat(HCr(:,:), Tr1(:,:), Tr2(:,:), Tr3(:,:), Tr4(:,:))
% [r,p]= corr(TdBu1,mRT,'rows','complete','type','spearman')
%%
% mRT1= mRT(40:69,:)
% mRT2= mRT(70:114,:)
% mRT3= mRT(115:148,:)
% mRT4= mRT(149:174,:)
% [r1,p1]= corr(Tr1,mRT1,'rows','complete','type','spearman')
% [r2,p2]= corr(Tr2,mRT2,'rows','complete','type','spearman')
% [r3,p3]= corr(Tr3,mRT3,'rows','complete','type','spearman')
% [r4,p4]= corr(Tr4,mRT4,'rows','complete','type','spearman')
% %%
% i= 1.5*iqr(HCr)
% Q1= quantile(HCr,.25)
% Q3= quantile(HCr,.75)
% F1= Q1-i
% F2= Q3+i
