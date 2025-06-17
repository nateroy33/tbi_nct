clear all
close all
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))

load mild_mRT_reorderd.mat
load correlation_metrics_mild_reordered.mat
condIND = [repelem(0,HCsubjs) repelem(1,T1subjs) repelem(2,T2subjs) repelem(3,T3subjs) repelem(4,T4subjs)];

MRT= mild_mRT_adjusted_z';

%%
%Correlation b/w mean TE (overall) and MRT
[HC_mean_z, mu, sigma] = zscore(avg_E(condIND==0));

TBI_mean_z = (avg_E(condIND~=0)-mu) / sigma;

avg_E_z2 = [HC_mean_z; TBI_mean_z];
avg_E = avg_E_z2;
[r3,p3]= corr(avg_E, MRT, 'rows', 'complete','type','Spearman')
cvec = [repmat([0 0 0],HCsubjs,1); repmat([.3010 .7450 .9330],T1subjs,1); repmat([.4940 .1840 .5560],T2subjs,1); repmat([.8500 .3250 .0980],T3subjs,1); repmat([.6350 .0780 .1840],T4subjs,1)];
%HC (no data), T1 red, T2 blue, T3 green, T4 black
figure; scatter(avg_E,MRT,[],cvec); lsline
title('Average Transition Energy vs Mean Reaction Time on the Attention Network Test','FontSize',18);
text(3, -3,['R = ',num2str(r3)],'FontSize',28)
text(3, -3.5,['p = ',num2str(p3)],'FontSize',28)
xlabel('Average Transition Energy','FontSize',18,'fontWeight','bold');
ylabel('Mean Reaction Time on the Attention Network Test','FontSize',18,'fontWeight','bold');
TickSize = 16
ax = gca
ax.FontSize = TickSize
% set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
% "higher avg TE = slower reaction time (poorer performace)
