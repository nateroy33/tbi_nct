clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
load TS86.mat
%%
HCsubjs= 39;
T1subjs= 30;
T2subjs= 45;
T3subjs= 34;
T4subjs= 26;
numsubjs= 174;
TR= 171;

subj= [HCsubjs T1subjs T2subjs T3subjs T4subjs];
maxsubj= max(subj);
%%
A = TS86

for i = 0:(numsubjs-1)
    
    A_Temp = A;
    
    rows = [TR*i+1,TR*i+TR];
    
    A_Temp = A_Temp((rows(1):rows(2)),:); %Temporary  171 x 86 matrix for input
        
    for j=1:width(A_Temp)
        
       s = sample_entropy(3,(0.2*std(A_Temp(:,j))), A_Temp(:,j), 1);
        
        Output(i+1,j) = s(1,1);
    end
        
end

VP = mean(Output,2)
%%

load editTP.mat unnamed1 % matrix specifying TP of TBI subjects
load CogDemo_174sub.mat
diag = diag'
scanInd0 = diag(:,:); % for HC
scanInd1 = unnamed1(:,1);% for TBI TP1
scanInd2 = unnamed1(:,2);% for TBI TP2
scanInd3 = unnamed1(:,3);% for TBI TP3
scanInd4 = unnamed1(:,4);% for TBI TP4

entropHC= VP(scanInd0==0,:)
entropHC= entropHC([1:21,23:end]);
entropT1= VP(scanInd1==1,:)
entropT2= VP(scanInd2==1,:)
entropT3= VP(scanInd3==1,:)
entropT4= VP(scanInd4==1,:)
%%
i= 1.5*iqr(entropHC)
F1= quantile(entropHC,.25)-i
F2= quantile(entropHC,.75)+i
%%
entrop = nan(maxsubj,5);
entrop(1:length(entropHC),1)= entropHC(:,:);
entrop(1:length(entropT1),2)= entropT1(:,:);
entrop(1:length(entropT2),3)= entropT2(:,:);
entrop(1:length(entropT3),4)= entropT3(:,:);
entrop(1:length(entropT4),5)= entropT4(:,:);
entropy= [entropT1; entropT2; entropT3; entropT4];

%%
figure;
violinplot(entrop);
title('Global Entropy');
xticklabels({'hc'; '1w post-TBI'; '1m post-TBI'; '6m post-TBI'; '12m post-TBI'});
% xlabel('post-TBI')
ylabel('Entropy');
TickSize = 16
ax = gca
ax.FontSize = TickSize
axis([0 6 0.88 1.10])
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
%%
[h1,p1,ci,tstat1]= ttest2(entropHC,entropT1)
[h1,p2,ci,tstat1]= ttest2(entropHC,entropT2)
[h1,p3,ci,tstat1]= ttest2(entropHC,entropT3)
[h1,p4,ci,tstat1]= ttest2(entropHC,entropT4)
pval= [p1 p2 p3 p4]
upval = mafdr(pval,'BHFDR',1);
%%
% [h2,p2,ci,stats]= ttest2(entrop(:,1),entrop(:,3))
% %%
% [h3,p3]= ttest2(entrop(:,1),entrop(:,4))
% [h4,p4]= ttest2(entrop(:,1),entrop(:,5))
%%
load gTEmildF.mat %nusubjs X (numClusters^2) TE matrix
m= mean(E_full, 2)
energyHC= m(scanInd0==0,:)
energyT1= m(scanInd1==1,:)
energyT2= m(scanInd2==1,:)
energyT3= m(scanInd3==1,:)
energyT4= m(scanInd4==1,:)

energy = nan(maxsubj,5);
energy(1:length(energyHC),1)= energyHC(:,:);
energy(1:length(energyT1),2)= energyT1(:,:);
energy(1:length(energyT2),3)= energyT2(:,:);
energy(1:length(energyT3),4)= energyT3(:,:);
energy(1:length(energyT4),5)= energyT4(:,:);

figure;

scatter(entrop(:,1),energy(:,1),"MarkerEdgeColor","w", "MarkerFaceColor",[0 0.4470 0.7410]);
title('Entropy vs. Energy')
xlabel('entropy')
ylabel('transition energy')
[r1,p12]= corr(m, VP,'rows','complete','type','Spearman');
text(.94, 55,['R = ',num2str(r1)]);
text(.94, 52,['p = ',num2str(p12)]);

hold on;
scatter(entrop(:,2),energy(:,2),"MarkerEdgeColor","w", "MarkerFaceColor",[0.3010 0.7450 0.9330]);

hold on;
scatter(entrop(:,3),energy(:,3), "MarkerEdgeColor","w", "MarkerFaceColor",[0.4940 0.1840 0.5560]);

hold on;
scatter(entrop(:,4),energy(:,4), "MarkerEdgeColor","w", "MarkerFaceColor",[0.8500 0.3250 0.0980]);

hold on;
scatter(entrop(:,5),energy(:,5), "MarkerEdgeColor","w", "MarkerFaceColor",[0.6350 0.0780 0.1840]);
%%
load mild_mRT_reorderd.mat
mRT=mRT'
mRT1= mRT(40:69,:)
 mRT2= mRT(70:114,:)
 mRT3= mRT(115:148,:)
 mRT4= mRT(149:174,:)
 [r1,p1]= corr(entropT1,mRT1,'rows','complete','type','spearman')
 [r2,p2]= corr(entropT2,mRT2,'rows','complete','type','spearman')
 [r3,p3]= corr(entropT3,mRT3,'rows','complete','type','spearman')
 [r4,p4]= corr(entropT4,mRT4,'rows','complete','type','spearman')
