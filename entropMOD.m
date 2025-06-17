clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
load TTS86.mat
%%
HC1subjs= 15;
HC2subjs= 11;
T1subjs= 42;
T2subjs= 17;
numsubjs= 85
TR= 420
subj= [HC1subjs HC2subjs T1subjs T2subjs];
maxsubj= max(subj);
%%
A = TTS86

for i = 0:(numsubjs-1)
    
    A_Temp = A;
    
    rows = [TR*i+1,TR*i+TR];
    
    A_Temp = A_Temp((rows(1):rows(2)),:); %Temporary  411 x 86 matrix for input
        
    for j=1:width(A_Temp)
        
       s = sample_entropy(3,(0.2*std(A_Temp(:,j))), A_Temp(:,j), 1);
        
        Output(i+1,j) = s(1,1);
    end
        
end

VP = mean(Output,2)

entropHC= VP(1:15,:)
entropT1= VP(27:68,:)
entropT2= VP(69:85,:)

i= 1.5*iqr(entropHC)
F1= quantile(entropHC,.25)-i
F2= quantile(entropHC,.75)+i

entrop = nan(42,3);
entrop(1:length(entropHC),1)= entropHC(:,:);
entrop(1:length(entropT1),2)= entropT1(:,:);
entrop(1:length(entropT2),3)= entropT2(:,:);

%%
figure;
violinplot(entrop);
title('Global Entropy');
xticklabels({'hc'; '4-6m post-TBI'; '12m post-TBI';});
% xlabel('post-TBI')
ylabel('Entropy');
TickSize = 16
ax = gca
ax.FontSize = TickSize
set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')
%xlim([0 4])
%ylim([0 0.6])

[h1,p1,ci,stat1]= ttest(entrop(:,1),entrop(:,2))
[h2,p2,ci,stat2]= ttest(entrop(:,1),entrop(:,3))
[h3,p3,ci,stat3]= ttest(entrop(:,2),entrop(:,3))
%%
% load E_fullF.mat
% m= mean(E_full, 2)
% energyHC= m(1:15,:)
% energyT1= m(27:67,:)
% energyT2= m(68:84,:)
% 
% energy = nan(41,3);
% energy(1:length(energyHC),1)= energyHC(:,:);
% energy(1:length(energyT1),2)= energyT1(:,:);
% energy(1:length(energyT2),3)= energyT2(:,:);

% figure;
% scatter(entrop(:,1),energy(:,1),"MarkerEdgeColor","w", "MarkerFaceColor",[0 0.4470 0.7410]);
% title('Entropy vs. Energy')
% xlabel('entropy')
% ylabel('transition energy')
% [r1,p12]= corr(m, VP,'rows','complete','type','Spearman');
% text(1, 140000,['R = ',num2str(r1)]);
% text(1, 120000,['p = ',num2str(p12)]);
% hold on;
% scatter(entrop(:,2),energy(:,2),"MarkerEdgeColor","w", "MarkerFaceColor",[0.3010 0.7450 0.9330]);
% hold on;
% scatter(entrop(:,3),energy(:,3), "MarkerEdgeColor","w", "MarkerFaceColor",[0.4940 0.1840 0.5560]);
