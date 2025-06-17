%% load stuff
clear all
close all

basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))

split = 3;
numClusters = 4;
TR = 420; 

load TTS86.mat
% load FODTAR2.4.mat 
TTS86= TTS86(1:35280,:); %LAST subject has no SC data
rTTS86= reshape(TTS86,420,84,86);
rTTS86= permute(rTTS86,[2,1,3]);
concTS = rTTS86;

nsubjs = size(rTTS86,1); 
HCsubjs= 15;
nparc = 86;
numclusters= 4;

subjInd=repelem(1:nsubjs,TR)';

clusterNames={};
clusterNamesUp={};
clusterNamesDown={};

%% compute centroids

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'partition','centroids'); 
partition= reshape(partition,[420,85])';
overallCentroids=centroids;

centroids=NaN(nsubjs,nparc,numClusters);

for i=1:nsubjs
    
    disp(['Generating centroids for scan ',num2str(i)]);
    centroids(i,:,:) = GET_CENTROIDS(squeeze(concTS(i,:,:)),partition(i,:),numClusters);
end
%% Top Down vs. Bottom Up
% load TdBu.mat %vector specifying top down vs. bottom up information for each ROI
% mcentroids= mean(centroids,1)
% mcentroids= squeeze(mcentroids)
% mcentroids= L2MAGNITUDENORM(mcentroids)
% mcentroids= abs(mcentroids)
% VisMD= dot(mcentroids(:,1), TdBu)
% DmnPD= dot(mcentroids(:,2), TdBu)
% SubMD= dot(mcentroids(:,3), TdBu)
% VisPD= dot(mcentroids(:,4), TdBu)
% %%
% c1= centroids(:,:,1); %reorder centroids based on top-down vs. bottom-up
% c2= centroids(:,:,2);
% c3= centroids(:,:,3);
% c4= centroids(:,:,4);
% centroids= cat(3,c2,c3,c4,c1)% 
%% load sc etc for E_calcs

load r3SCmatsTBI.mat
SCmatsZ= SCmats(:,:,1);
SCmatsZ(find(eye(nparc,nparc)==1))=0;
for i= 2:84
SCmat= SCmats(:,:,i);
SCmat(find(eye(nparc,nparc)==1))=0;
SCmatsZ= cat(3,SCmatsZ,SCmat);
end
SCmats= SCmatsZ;
a= SCmats(:,:,1)
HCSCmats= SCmats(:,:,1:15);
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

E_full=NaN(nsubjs,numClusters^2);

for i=1:nsubjs
    sci = squeeze(SCmats(:,:,i));
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

[h1,p1,ci1,tstat1] = ttest2(E_full(1:15,:),E_full(27:67,:));
[h2,p2,ci2,tstat2] = ttest2(E_full(1:15,:),E_full(68:84,:));
[h3,p3,ci3,tstat3] = ttest2(E_full(27:67,:),E_full(68:84,:));

figure;
imagesc([tstat1.tstat; tstat2.tstat; tstat3.tstat])
allt = [tstat1.tstat; tstat2.tstat; tstat3.tstat];
pall = [p1; p2; p3];
% p1= reshape(p1,[numClusters numClusters])';
p1u= mafdr(p1,'BHFDR',1)
p1= reshape(p1,[numClusters numClusters])';
p1u= reshape(p1u,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
p2u= mafdr(p2,'BHFDR',1)
p2= reshape(p2,[numClusters numClusters])';
p2u= reshape(p2u,[numClusters numClusters])';
% p3= reshape(p2,[numClusters numClusters])';
p3u= mafdr(p3,'BHFDR',1)

maxVal = max(allt(:)); % sync color scales
minVal = min(allt(:));
sig_thresh = 0.05;
clusterNames = {'DMN+','SUB-','VIS+','VIS-'};

subplot(1,4,1);
imagesc(reshape(tstat1.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 4-6 months post-TBI');
set(gca,'FontSize',12);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar
hold on; 
sig_thresh = 0.05; 
[y,x] = find(squeeze(p1u) < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(p1) < sig_thresh); 
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
% [y,x] = find(p1 < sig_thresh); 
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% [y2,x2] = find(squeeze(p1u) < sig_thresh); 
% if ~isempty(y2)
%     text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
% end
subplot(1,4,2);
imagesc(reshape(tstat2.tstat,[numClusters numClusters])',[minVal maxVal])
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('HC vs 12 months post-TBI');
set(gca,'FontSize',12);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar
sig_thresh = 0.05; 
[y,x] = find(squeeze(p2u) < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 20);
[y2,x2] = find(squeeze(p2) < sig_thresh); 
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 20);
% [y,x] = find(p2 < sig_thresh); 
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
% [y2,x2] = find(p2u < sig_thresh); 
% if ~isempty(y2)
%     text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
% end

% subplot(1,4,3);
% imagesc(reshape(tstat3.tstat,[numClusters numClusters])',[minVal maxVal])
% xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('Initial State'); xlabel('Final State');
% title('4-6 months vs 12 months post-TBI');
% set(gca,'FontSize',12);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
% caxis([minVal maxVal]); colorbar
% [y,x] = find(p2 < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
% if ~isempty(y)
%     text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
% end
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
% 
% [h1,p1,ci1,tstat1] = ttest2(E_full(1:11,:),E_full(22:49,:));
% [h2,p2,ci2,tstat2] = ttest2(E_full(1:11,:),E_full(50:60,:));
% 
% figure;
% imagesc([tstat1.tstat; tstat2.tstat])
% allt = [tstat1.tstat; tstat2.tstat];
% pall = [p1; p2];
% p1= reshape(p1,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
% 
% maxVal = max(allt(:)); % sync color scales
% minVal = min(allt(:));
% sig_thresh = 0.05;
% clusterNames = {'VIS+','VIS-','DMN+','SOM-'};
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
% %% Flavor 3: Indiv SC/Grp Avg Centroids
% E_full=NaN(nsubjs,numClusters^2);
% 
% for i=1:nsubjs
%     sci = squeeze(SCmats(:,:,i));
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
% 
% [h1,p1,ci1,tstat1] = ttest2(E_full(1:11,:),E_full(22:49,:));
% [h2,p2,ci2,tstat2] = ttest2(E_full(1:11,:),E_full(50:60,:));
% 
% figure;
% imagesc([tstat1.tstat; tstat2.tstat])
% allt = [tstat1.tstat; tstat2.tstat];
% pall = [p1; p2];
% p1= reshape(p1,[numClusters numClusters])';
% p2= reshape(p2,[numClusters numClusters])';
% 
% maxVal = max(allt(:)); % sync color scales
% minVal = min(allt(:));
% sig_thresh = 0.05;
% clusterNames = {'VIS+','VIS-','DMN+','SOM-'};
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

%% Avg E by Time Point

[h1,p1,ci1,tstat1] = ttest2(nanmean(E_full(1:15,:),2),nanmean(E_full(27:67,:),2));
[h2,p2,ci2,tstat2] = ttest2(nanmean(E_full(1:15,:),2),nanmean(E_full(68:84,:),2));
pval=[p1 p2]
pvalu= mafdr(pval,'BHFDR',1)

AEhc= nanmean(E_full(1:15,:),2)
i= 1.5*iqr(AEhc)
Q1= quantile(AEhc,.25)
Q3= quantile(AEhc,.75)
F1= Q1-i
F2= Q3+i

data= NaN(41,3);
data(1:15,1)= nanmean(E_full(1:15,:),2); %HC
data(1:41,2)= nanmean(E_full(27:67,:),2); %T1
data(1:17,3)= nanmean(E_full(68:84,:),2); %T2
AEmod= cat(1,(nanmean(E_full(1:15,:),2)),(nanmean(E_full(27:67,:),2)),(nanmean(E_full(68:84,:),2)))

figure;
violinplot(data);
title('Average Transition Energy by Group');
ylabel('Average Transition Energy');
xticklabels({'hc'; '4-6m post-TBI'; '12m post-TBI'});
%xlabel('post-TBI')
TickSize = 16
ax = gca
ax.FontSize = TickSize
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

[h,p1,Ci,stats]= ttest2(data(:,1),data(:,2))
[h,p2,Ci,stats]= ttest2(data(:,1),data(:,3))
pval = [p1 p2]
upval = mafdr(pval,'BHFDR',1)
% %% By Column
% 
% VisM= E_full(:,[2 6 10 14]);
% VisP= E_full(:,[1 5 9 13]);
% DmnP= E_full(:,[3 7 11 15]);
% SomM= E_full(:,[4 8 12 16]);
% 
% 
% % [h1,p1,ci1,tstat1] = ttest2(nanmean(DmnP(1:11,:),2),nanmean(DmnP(22:49,:),2));
% % [h2,p2,ci2,tstat2] = ttest2(nanmean(DmnP(1:11,:),2),nanmean(DmnP(50:60,:),2));
% 
% data= NaN(28,3);
% data(1:11,1)= nanmean(E_full(1:11,:),2);
% data(1:28,2)= nanmean(E_full(22:49,:),2);
% data(1:11,3)= nanmean(E_full(50:60,:),2);
% 
% figure;
% violinplot(data);
% title('Avg E for Transitions Into DMN+ by Group');
% ylabel('Avg E');
% xticklabels({'hc'; 'tbi1'; 'tbi2'});


% % New order (VisM, DmnP, SubM, VisP)
% 
% E_hc= E_full(1:15,:)
% E_hcm= (E_hc)
% HCbu= (E_hcm(:,2)+E_hcm(:,3)+E_hcm(:,4)+E_hcm(:,7)+E_hcm(:,8)+E_hcm(:,12))/6 %(trans into more bu state)
% HCtd= (E_hcm(:,5)+E_hcm(:,9)+E_hcm(:,10)+E_hcm(:,13)+E_hcm(:,14)+E_hcm(:,15))/6 %(trans into more td state)
% [h,p,Ci1,stats]= ttest(HCbu,HCtd)
% %%
% E_1= E_full(27:67,:)
% E_1m= (E_1)
% T1bu= (E_1m(:,2)+E_1m(:,3)+E_1m(:,4)+E_1m(:,7)+E_1m(:,8)+E_1m(:,12))/6 %(trans into more bu state)
% T1td= (E_1m(:,5)+E_1m(:,9)+E_1m(:,10)+E_1m(:,13)+E_1m(:,14)+E_1m(:,15))/6 %(trans into more td state)
% [h,p,Ci2,stats]= ttest(T1bu,T1td)
% %%
% E_2= E_full(68:84,:)
% E_2m= (E_2)
% T2bu= (E_2m(:,2)+E_2m(:,3)+E_2m(:,4)+E_2m(:,7)+E_2m(:,8)+E_2m(:,12))/6 %(trans into more bu state)
% T2td= (E_2m(:,5)+E_2m(:,9)+E_2m(:,10)+E_2m(:,13)+E_2m(:,14)+E_2m(:,15))/6 %(trans into more td state)
% [h,p,Ci3,stats]= ttest(T2bu,T2td)
% %%
% % HCbu= (E_hcm(:,2)+E_hcm(:,3)+E_hcm(:,4)+E_hcm(:,7)+E_hcm(:,8)+E_hcm(:,12)) %(trans into more bu state)
% % HCtd= (E_hcm(:,5)+E_hcm(:,9)+E_hcm(:,10)+E_hcm(:,13)+E_hcm(:,14)+E_hcm(:,15)) %(trans into more td state)
% HCr= HCbu./HCtd
% 
% % T1bu= (E_1m(:,2)+E_1m(:,3)+E_1m(:,4)+E_1m(:,7)+E_1m(:,8)+E_1m(:,12)) %(trans into more bu state)
% % T1td= (E_1m(:,5)+E_1m(:,9)+E_1m(:,10)+E_1m(:,13)+E_1m(:,14)+E_1m(:,15)) %(trans into more td state)
% Tr1= T1bu./T1td
% 
% % T2bu= (E_2m(:,2)+E_2m(:,3)+E_2m(:,4)+E_2m(:,7)+E_2m(:,8)+E_2m(:,12)) %(trans into more bu state)
% % T2td= (E_2m(:,5)+E_2m(:,9)+E_2m(:,10)+E_2m(:,13)+E_2m(:,14)+E_2m(:,15)) %(trans into more td state)
% Tr2= T2bu./T2td
% 
% data33= NaN(41,3);
% data33(1:15,1)= HCr;
% data33(1:41,2)= Tr1;
% data33(1:17,3)= Tr2;
% 
% 
% figure;
% violinplot(data33);
% title('Transition energy ratio: Top Down vs. Bottom Up');
% ylabel('Top Down Energy/Bottom Up Energy');
% xlabel('post-TBI');
% xticklabels({'hc'; '4-6 months'; '12 months'});
% ax = gca
% ax.FontSize = TickSize
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
% % text(4, 1.5,['*'],'FontSize',28)
% % text(5, 1.5,['*'], 'FontSize',28)
% 
% [h,p,Ci,stats]= ttest2(data33(:,2),data33(:,3))
% 
%% New TdBu
E_hc= E_full(1:15,:);
E_1= E_full(27:67,:);
E_2= E_full(68:84,:);

HCr= (((E_hc(:,5)-E_hc(:,2))+(E_hc(:,10)-E_hc(:,7))+(E_hc(:,15)-E_hc(:,12))+(E_hc(:,9)-E_hc(:,3))+(E_hc(:,14)-E_hc(:,8))+(E_hc(:,13)-E_hc(:,4)))./6)
Tr1= (((E_1(:,5)-E_1(:,2))+(E_1(:,10)-E_1(:,7))+(E_1(:,15)-E_1(:,12))+(E_1(:,9)-E_1(:,3))+(E_1(:,14)-E_1(:,8))+(E_1(:,13)-E_1(:,4)))./6)
Tr2= (((E_2(:,5)-E_2(:,2))+(E_2(:,10)-E_2(:,7))+(E_2(:,15)-E_2(:,12))+(E_2(:,9)-E_2(:,3))+(E_2(:,14)-E_2(:,8))+(E_2(:,13)-E_2(:,4)))./6)

data33= NaN(41,3);
data33(1:15,1)= HCr;
data33(1:41,2)= Tr1;
data33(1:17,3)= Tr2;

figure;
violinplot(data33);
title('Transition Energy Asymmetry');
ylabel('Bottom Up Energy - Top Down Energy');
%xlabel('post-TBI');
TickSize = 16
ax = gca
ax.FontSize = TickSize
xticklabels({'hc'; '4-6m post-TBI'; '12m post-TBI'});
text(4, 1.5,['*'],'FontSize',28)
text(5, 1.5,['*'], 'FontSize',28)
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')

%%
[h,p,Ci,stats]= ttest2(HCr,Tr1)