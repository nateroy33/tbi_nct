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



%% New Top Down / Bottom Up E Calc
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

