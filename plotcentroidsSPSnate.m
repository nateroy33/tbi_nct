%% plot centroids on gummibrain from Keith Jamison-> https://github.com/kjamison/atlasblobs

clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);	


%% plot 

numClusters = 6;
split= 1

load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

load(fullfile('atlasblobs_saved.mat'))

%choose atlas
%options for whichatlas: aal, cc200, cc400, ez, ho, tt, fs86, ls463,
%sch232, sch454
whichatlas={'fs86'};

%%

%set colormap
cmap=plasma;


%%
%set data you want to plot
% load m.mat; %average centroids (nparc X numclusters)
% data = centroids;
% data = zscore(mean_node_diff');
% Acentroids= squeeze(Acentroids);
data = centroids;
%set min/max limits for plot
% clim=[-1 1];
clim= [min(data,[],'all') max(data,[],'all')]

cmap = plasma;
% 
f=figure;


for i=1:numClusters
    
    subplot(1,numClusters,i);
    
    img=display_atlas_blobs(data(:,i),atlasblobs,...
        'atlasname',whichatlas,...
        'render',true,...
        'backgroundimage',true,...
        'crop',true,...
        'colormap',cmap,...
        'clim', clim,...
        'backgroundcolor',[1 1 1]);%,...
%         'roimask',roimask); %last argmument optional, depends if you need to remove roi's
    
    
    imshow(img);
    c=colorbar('SouthOutside', 'fontsize', 16);
%     c.Label.String='Low-Amplitude Activity  High Amplitude Activity';
    set(gca,'colormap',cmap);
    caxis(clim);
%     title(clusterNames{i});
    
%     annotation('textbox',[.89 .33 .1 .2],'String','RH','EdgeColor','none','fontname','Arial','fontsize',16,'color','black')
%     annotation('textbox',[.07 .33 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','black')
%     annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none','fontname','Arial', 'fontsize',16,'color','black', 'horizontalalignment','center')
%     annotation('textbox',[.45 .1 .1 .04],'String','Medial','EdgeColor','none','fontname','Arial', 'fontsize',16,'color','black', 'horizontalalignment','center')
%     
end

% figure;
% title(clusterNames{i});
% 