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


