% a=clock;
% rng(a(6));

% load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
E = load('TTS86.mat');
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
savedir = fullfile(basedir,'repkmeans2');
mkdir(savedir); cd(savedir);
concTS = E.TTS86; % csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
nreps= 50; 
maxI= 1000;
distanceMethod = 'correlation';
%%% MAKE SURE TO RUN THE CORRECT SPLIT
split= 3; 
% split 1 is the data as received from Amy
% split 2 is data received from Keith in 10-2021

% if zdim > 0
% 	concTS = zscore(concTS,[],zdim);
% end

disp('data loaded');
for numClusters= 2:13
    
    disp(['K = ',num2str(numClusters),'Split = ',num2str(split)]);
    
%     savedir = ['kmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
%     mkdir(savedir);
%     %%
%     cd(savedir);
    
    disp('start k-means');
%     for R = 1:nreps
        [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates', nreps, 'MaxIter', maxI);
        save(['kmeans',num2str(split),'k_',num2str(numClusters),'.mat'],'partition','sumd');
        clear partition
        clear sumd
%        disp(['K-means ',num2str(R)]);
%     end
end
disp('complete');
