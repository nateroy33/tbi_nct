clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))

load Violin5Data_bp1_k4.mat
HCsubjs= 39;
T1subjs= 30;
T2subjs= 45;
T3subjs= 34;
T4subjs= 26;

subj= [HCsubjs T1subjs T2subjs T3subjs T4subjs];
maxsubj= max(subj);
%Fractional Occupancy
HCfo= HCfo';
T1fo= T1fo';
T2fo= T2fo';
T3fo= T3fo';
T4fo= T4fo';

FO1= nan(maxsubj,5);
FO1(1:length(HCfo),1)= HCfo(:,1);
FO1(1:length(T1fo),2)= T1fo(:,1);
FO1(1:length(T2fo),3)= T2fo(:,1);
FO1(1:length(T3fo),4)= T3fo(:,1);
FO1(1:length(T4fo),5)= T4fo(:,1);
%%
a= anova1(FO1);
[h,p1,ci,stats] = ttest2(FO1(:,1),FO1(:,2));
%%
[h,p2,ci,stats] = ttest2(FO1(:,1),FO1(:,3));
%%
[h,p3,ci,stats] = ttest2(FO1(:,1),FO1(:,4));
%%
[h,p4,ci,stats] = ttest2(FO1(:,1),FO1(:,5));
%%
FO2= nan(maxsubj,5);
FO2(1:length(HCfo),1)= HCfo(:,2);
FO2(1:length(T1fo),2)= T1fo(:,2);
FO2(1:length(T2fo),3)= T2fo(:,2);
FO2(1:length(T3fo),4)= T3fo(:,2);
FO2(1:length(T4fo),5)= T4fo(:,2);
%%
c= anova1(FO2);
[h,p5,ci,stats] = ttest2(FO2(:,1),FO2(:,2));
%%
[h,p6,ci,stats] = ttest2(FO2(:,1),FO2(:,3));
%%
[h,p7,ci,stats] = ttest2(FO2(:,1),FO2(:,4));
%%
[h,p8,ci,stats] = ttest2(FO2(:,1),FO2(:,5));
%%
b= anova1(FO2);

FO3= nan(maxsubj,5);
FO3(1:length(HCfo),1)= HCfo(:,3);
FO3(1:length(T1fo),2)= T1fo(:,3);
FO3(1:length(T2fo),3)= T2fo(:,3);
FO3(1:length(T3fo),4)= T3fo(:,3);
FO3(1:length(T4fo),5)= T4fo(:,3);
%%
c= anova1(FO3);
[h,p9,ci,stats] = ttest2(FO3(:,1),FO3(:,2));
%%
[h,p10,ci,stats] = ttest2(FO3(:,1),FO3(:,3));
%%
[h,p11,ci,stats] = ttest2(FO3(:,1),FO3(:,4));
%%
[h,p12,ci,stats] = ttest2(FO3(:,1),FO3(:,5));
%%
FO4= nan(maxsubj,5);
FO4(1:length(HCfo),1)= HCfo(:,4);
FO4(1:length(T1fo),2)= T1fo(:,4);
FO4(1:length(T2fo),3)= T2fo(:,4);
FO4(1:length(T3fo),4)= T3fo(:,4);
FO4(1:length(T4fo),5)= T4fo(:,4);
%%
d= anova1(FO4);
[h,p13,ci,stats] = ttest2(FO4(:,1),FO4(:,2));
%%
[h,p14,ci,stats] = ttest2(FO4(:,1),FO4(:,3));
%%
[h,p15,ci,stats] = ttest2(FO4(:,1),FO4(:,4));
%%
[h,p16,ci,stats] = ttest2(FO4(:,1),FO4(:,5));
%
% FO5= nan(maxsubj,5);
% FO5(1:length(HCfo),1)= HCfo(:,5);
% FO5(1:length(T1fo),2)= T1fo(:,5);
% FO5(1:length(T2fo),3)= T2fo(:,5);
% FO5(1:length(T3fo),4)= T3fo(:,5);
% FO5(1:length(T4fo),5)= T4fo(:,5);
%%
% d= anova1(FO5);
% [h,p17,ci,stats] = ttest2(FO5(:,1),FO5(:,2));
% %%
% [h,p18,ci,stats] = ttest2(FO5(:,1),FO5(:,3));
% %%
% [h,p19,ci,stats] = ttest2(FO5(:,1),FO5(:,4));
% %%
% [h,p20,ci,stats] = ttest2(FO5(:,1),FO5(:,5));
% %%
% FO6= nan(maxsubj,6);
% FO6(1:length(HCfo),1)= HCfo(:,6);
% FO6(1:length(T1fo),2)= T1fo(:,6);
% FO6(1:length(T2fo),3)= T2fo(:,6);
% FO6(1:length(T3fo),4)= T3fo(:,6);
% FO6(1:length(T4fo),5)= T4fo(:,6);
% %
% d= anova1(FO6);
% [h,p21,ci,stats] = ttest2(FO6(:,1),FO6(:,2));
% %%
% [h,p22,ci,stats] = ttest2(FO6(:,1),FO6(:,3));
% %%
% [h,p23,ci,stats] = ttest2(FO6(:,1),FO6(:,4));
% %%
% [h,p24,ci,stats] = ttest2(FO6(:,1),FO6(:,5));
% pval= [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24]
% upval = mafdr(pval,'BHFDR',1);
% comp = [pval; upval];
%%
figure;
subplot(1,4,1);
violinplot(FO1);
subtitle('DMN+/FPN-');
xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% xlabel('post-TBI')
TickSize = 16
ax = gca
ax.FontSize = TickSize
ylabel('Fractional Occupancy');
ylim([0 0.5])
subplot(1,4,2);
violinplot(FO2);
subtitle('DMN-/VATSOM+');
xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% xlabel('post-TBI')
ax = gca
ax.FontSize = TickSize
ylabel('Fractional Occupancy');
ylim([0 0.5])
subplot(1,4,3);
violinplot(FO3);
subtitle('DMN+/VATSOM-');
xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% xlabel('post-TBI')
ax = gca
ax.FontSize = TickSize
ylabel('Fractional Occupancy');
ylim([0 0.5])
subplot(1,4,4);
violinplot(FO4);
subtitle('DMN-/FPN+'); 
xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% xlabel('post-TBI')
ax = gca
ax.FontSize = TickSize
ylabel('Fractional Occupancy');
ylim([0 0.5])
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
% subplot(1,6,5);
% violinplot(FO5);
% subtitle('VIS+'); 
% xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% %xlabel('post-TBI')
% ax = gca
% ax.FontSize = TickSize
% ylabel('Fractional Occupancy');
% ylim([0 0.5])
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
% subplot(1,6,6);
% violinplot(FO6);
% subtitle('VIS-'); 
% xticklabels({'hc','1w post-TBI','1m post-TBI','6m post-TBI','12m post-TBI'});
% % xlabel('post-TBI')
% ax = gca
% ax.FontSize = TickSize
% ylabel('Fractional Occupancy');
% ylim([0 0.5])
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
%%
% %% Dwell Time
% HCdt= HCdt';
% T1dt= T1dt';
% T2dt= T2dt';
% T3dt= T3dt';
% T4dt= T4dt';
% 
% DT1= nan(maxsubj,5);
% DT1(1:length(HCdt),1)= HCdt(:,1);
% DT1(1:length(T1dt),2)= T1dt(:,1);
% DT1(1:length(T2dt),3)= T2dt(:,1);
% DT1(1:length(T3dt),4)= T3dt(:,1);
% DT1(1:length(T4dt),5)= T4dt(:,1);
% 
% e= anova1(DT1);
% 
% DT2= nan(maxsubj,5);
% DT2(1:length(HCdt),1)= HCdt(:,2);
% DT2(1:length(T1dt),2)= T1dt(:,2);
% DT2(1:length(T2dt),3)= T2dt(:,2);
% DT2(1:length(T3dt),4)= T3dt(:,2);
% DT2(1:length(T4dt),5)= T4dt(:,2);
% 
% f= anova1(DT2);
% 
% DT3= nan(maxsubj,5);
% DT3(1:length(HCdt),1)= HCdt(:,3);
% DT3(1:length(T1dt),2)= T1dt(:,3);
% DT3(1:length(T2dt),3)= T2dt(:,3);
% DT3(1:length(T3dt),4)= T3dt(:,3);
% DT3(1:length(T4dt),5)= T4dt(:,3);
% 
% g= anova1(DT3);
% [h,p13,ci,stats] = ttest2(DT3(:,1),DT3(:,2));
% [h,p14,ci,stats] = ttest2(DT3(:,1),DT3(:,3));
% %%
% [h,p15,ci,stats] = ttest2(DT3(:,1),DT3(:,4));
% %%
% [h,p16,ci,stats] = ttest2(DT3(:,1),DT3(:,5));
% 
% DT4= nan(maxsubj,5);
% DT4(1:length(HCdt),1)= HCdt(:,4);
% DT4(1:length(T1dt),2)= T1dt(:,4);
% DT4(1:length(T2dt),3)= T2dt(:,4);
% DT4(1:length(T3dt),4)= T3dt(:,4);
% DT4(1:length(T4dt),5)= T4dt(:,4);
% 
% hh= anova1(DT4);

% figure;
% 
% subplot(1,4,1);
% violinplot(DT1);
% subtitle('MS-1a');
% subplot(1,4,2);
% violinplot(DT2);
% subtitle('MS-1b');
% subplot(1,4,3);
% violinplot(DT3);
% subtitle('MS-2a *');
% subplot(1,4,4);
% violinplot(DT4);
% subtitle('MS-2b');
% title('Dwell Time')
% 
% %Appearance Rate
% HCar= HCar';
% T1ar= T1ar';
% T2ar= T2ar';
% T3ar= T3ar';
% T4ar= T4ar';
% 
% AR1= nan(maxsubj,5);
% AR1(1:length(HCar),1)= HCar(:,1);
% AR1(1:length(T1ar),2)= T1ar(:,1);
% AR1(1:length(T2ar),3)= T2ar(:,1);
% AR1(1:length(T3ar),4)= T3ar(:,1);
% AR1(1:length(T4ar),5)= T4ar(:,1);
% 
% i= anova1(AR1);
% [h,p17,ci,stats] = ttest2(AR1(:,1),AR1(:,2));
% [h,p18,ci,stats] = ttest2(AR1(:,1),AR1(:,3));
% [h,p19,ci,stats] = ttest2(AR1(:,1),AR1(:,4));
% [h,p20,ci,stats] = ttest2(AR1(:,1),AR1(:,5));
% 
% AR2= nan(maxsubj,5);
% AR2(1:length(HCar),1)= HCar(:,2);
% AR2(1:length(T1ar),2)= T1ar(:,2);
% AR2(1:length(T2ar),3)= T2ar(:,2);
% AR2(1:length(T3ar),4)= T3ar(:,2);
% AR2(1:length(T4ar),5)= T4ar(:,2);
% 
% j= anova1(AR2);
% 
% AR3= nan(maxsubj,5);
% AR3(1:length(HCar),1)= HCar(:,3);
% AR3(1:length(T1ar),2)= T1ar(:,3);
% AR3(1:length(T2ar),3)= T2ar(:,3);
% AR3(1:length(T3ar),4)= T3ar(:,3);
% AR3(1:length(T4ar),5)= T4ar(:,3);
% 
% k= anova1(AR3);
% 
% AR4= nan(maxsubj,5);
% AR4(1:length(HCar),1)= HCar(:,4);
% AR4(1:length(T1ar),2)= T1ar(:,4);
% AR4(1:length(T2ar),3)= T2ar(:,4);
% AR4(1:length(T3ar),4)= T3ar(:,4);
% AR4(1:length(T4ar),5)= T4ar(:,4);
% 
% l= anova1(AR4);
% 
% figure;
% subplot(1,4,1);
% violinplot(AR1);
% subtitle('MS-1a *');
% subplot(1,4,2);
% violinplot(AR2);
% subtitle('MS-1b');
% subplot(1,4,3);
% violinplot(AR3);
% subtitle('MS-2a');
% subplot(1,4,4);
% violinplot(AR4);
% subtitle('MS-2b');
% title('Appearance Rate');
