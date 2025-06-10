clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 6;
split=3; 
savedir = fullfile(basedir,'results','example2');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
TR = 420;
[X,~]=size(partition);

%% Reassign clusters in partition 
%%
% clusterNames = {'MS-1';'MS-2'}; % new combined cluster names
clusterNames = {'LIM+/FPN-';'DMN-/FPN-';'DMN+/FPN+';'LIM-/FPN+';'VIS+';'VIS-'}%;'MS-3'}%'MS-4a';'MS-4b'};
% numClusters = 4; %new number of clusters (combined)
%%
c_partition = NaN(X,1); %initialize new partition (combined part.)

c_partition(partition == 1) = 5; 
c_partition(partition == 2) = 2; 

c_partition(partition == 3) = 3; 
c_partition(partition == 4) = 4; %reassign all of # in parentheses to be the number on the right

 c_partition(partition == 5) = 6; 
 c_partition(partition == 6) = 1;
partition = c_partition; 

%% reorder centroids
[nparc,~]=size(centroids);
r_centroids = NaN(nparc,numClusters);

r_centroids(:,1) = centroids(:,6); % new (# on left) is old (# on right) FLIPPED ORDER FROM ABOVE SECTION
r_centroids(:,2) = centroids(:,2);
r_centroids(:,3) = centroids(:,3);
r_centroids(:,4) = centroids(:,4);
r_centroids(:,5) = centroids(:,1);
r_centroids(:,6) = centroids(:,5);

centroids = r_centroids;

%% save partition in separate file


save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids','D');
