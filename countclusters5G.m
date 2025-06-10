clear all; close all;clc
basedir = 'C:\Users\nroy2\OneDrive - Cornell University\Documents\MATLAB\New\brain_states-master';
cd(basedir);
addpath(genpath('code'))
load tbi_cat.mat

%% Set inputs
split=1; %Split 1 is mild TBI set
numClusters= 6;
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
TR = 171;
load CogDemo_174sub.mat diag % Vector specifiying HC vs. TBI
load editTP.mat unnamed1 % Vector specifying timepoint of each TBI subject
TPtbi= unnamed1;
scanInd0 = repelem(diag,TR); % for HC partition
scanInd1 = repelem(TPtbi(:,1),TR); % for TBI1 partition
scanInd2 = repelem(TPtbi(:,2),TR); % for TBI2 partition
scanInd3 = repelem(TPtbi(:,3),TR); % for TBI3 partition
scanInd4 = repelem(TPtbi(:,4),TR); % for TBI4 partition
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

% Partitions of time series
HCPartition= partition(scanInd0 == 0);
T1Partition= partition(scanInd1 == 1);
T2Partition= partition(scanInd2 == 1);
T3Partition= partition(scanInd3 == 1);
T4Partition= partition(scanInd4 == 1);

HCsubjs= 39;
HCsubjInd= repelem(1:HCsubjs, TR)';

T1subjs = 30;
T1subjInd= repelem(1:T1subjs, TR)';

T2subjs = 45;
T2subjInd= repelem(1:T2subjs, TR)';

T3subjs = 34;
T3subjInd= repelem(1:T3subjs, TR)';

T4subjs = 26;
T4subjInd= repelem(1:T4subjs, TR)';

scan_length = 7; %time of each scan (in min)
rep_time = 2; %length of each TR (in s)


%% big loop
for numClusters=[6]
    load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
    %% count number of each cluster per scan
    
    HCP=reshape(HCPartition,TR,HCsubjs);
    HCcount = zeros(numClusters,HCsubjs);
    for b=1:HCsubjs
        [HCcount(:,b),~] = hist(HCP(:,b),1:numClusters);
    end
    
    T1P=reshape(T1Partition,TR,T1subjs);
    T1count = zeros(numClusters,T1subjs);
    for b=1:T1subjs
        [T1count(:,b),~] = hist(T1P(:,b),1:numClusters);
    end
    
    T2P=reshape(T2Partition,TR,T2subjs);
    T2count = zeros(numClusters,T2subjs);
    for b=1:T2subjs
        [T2count(:,b),~] = hist(T2P(:,b),1:numClusters);
    end
    
    T3P=reshape(T3Partition,TR,T3subjs);
    T3count = zeros(numClusters,T3subjs);
    for b=1:T3subjs
        [T3count(:,b),~] = hist(T3P(:,b),1:numClusters);
    end
    
    T4P=reshape(T4Partition,TR,T4subjs);
    T4count = zeros(numClusters,T4subjs);
    for b=1:T4subjs
        [T4count(:,b),~] = hist(T4P(:,b),1:numClusters);
    end
    %% Calculate Fractional Occupancy
    HCfo = HCcount/TR; %data for all scans
    T1fo = T1count/TR;
    T2fo = T2count/TR;
    T3fo = T3count/TR;
    T4fo = T4count/TR;
    

%% Calculate Dwell Time and Appearance Rate
    %%HC
    HCdwell = []; 
    sumApp=0;
    for c=1:numClusters
        for b=1:HCsubjs
            appear=find(HCP(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    HCdwell(c,i,b)=s;
                    a=a+1;
                else
                    HCdwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(HCdwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    HC_count=zeros(numClusters,3,HCsubjs); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    HCdt=[];
    HCar=[];
    [~,maxIndex,~]=size(HCdwell);
    for c=1:numClusters
        for b=1:HCsubjs
            HC_count(c,1,b)=HC_count(c,1,b)+sum(HCdwell(c,:,b));
            a=1;
            while a<=maxIndex
                if HCdwell(c,a,b)~= 0
                    HC_count(c,2,b)=HC_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            HCdt(c,b)=HC_count(c,1,b)/HC_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            HCar(c,b)=HC_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
    

    %T1
    T1dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:T1subjs
            appear=find(T1P(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    T1dwell(c,i,b)=s;
                    a=a+1;
                else
                    T1dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(T1dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    T1_count=zeros(numClusters,3,T1subjs); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    T1dt=[];
    T1ar=[];
    [~,maxIndex,~]=size(T1dwell);
    for c=1:numClusters
        for b=1:T1subjs
            T1_count(c,1,b)=T1_count(c,1,b)+sum(T1dwell(c,:,b));
            a=1;
            while a<=maxIndex
                if T1dwell(c,a,b)~= 0
                    T1_count(c,2,b)=T1_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            T1dt(c,b)=T1_count(c,1,b)/T1_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            T1ar(c,b)=T1_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
  

     T2dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:T2subjs
            appear=find(T2P(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    T2dwell(c,i,b)=s;
                    a=a+1;
                else
                    T2dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(T2dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    T2_count=zeros(numClusters,3,T2subjs); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    T2dt=[];
    T2ar=[];
    [~,maxIndex,~]=size(T2dwell);
    for c=1:numClusters
        for b=1:T2subjs
            T2_count(c,1,b)=T2_count(c,1,b)+sum(T2dwell(c,:,b));
            a=1;
            while a<=maxIndex
                if T2dwell(c,a,b)~= 0
                    T2_count(c,2,b)=T2_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            T2dt(c,b)=T2_count(c,1,b)/T2_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            T2ar(c,b)=T2_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
  
    
     T3dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:T3subjs
            appear=find(T3P(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    T3dwell(c,i,b)=s;
                    a=a+1;
                else
                    T3dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(T3dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    T3_count=zeros(numClusters,3,T3subjs); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    T3dt=[];
    T3ar=[];
    [~,maxIndex,~]=size(T3dwell);
    for c=1:numClusters
        for b=1:T3subjs
            T3_count(c,1,b)=T3_count(c,1,b)+sum(T3dwell(c,:,b));
            a=1;
            while a<=maxIndex
                if T3dwell(c,a,b)~= 0
                    T3_count(c,2,b)=T3_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            T3dt(c,b)=T3_count(c,1,b)/T3_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            T3ar(c,b)=T3_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
  
    
     T4dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:T4subjs
            appear=find(T4P(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    T4dwell(c,i,b)=s;
                    a=a+1;
                else
                    T4dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(T4dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    T4_count=zeros(numClusters,3,T4subjs); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    T4dt=[];
    T4ar=[];
    [~,maxIndex,~]=size(T4dwell);
    for c=1:numClusters
        for b=1:T4subjs
            T4_count(c,1,b)=T4_count(c,1,b)+sum(T4dwell(c,:,b));
            a=1;
            while a<=maxIndex
                if T4dwell(c,a,b)~= 0
                    T4_count(c,2,b)=T4_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            T4dt(c,b)=T4_count(c,1,b)/T4_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            T4ar(c,b)=T4_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
  

    %% compute DT's and AR's averaged over clusters (does LSD overall have
    %lower DT's and higher AR's?)
    
    hcDT=mean(HCdt,1);
    t1DT=mean(T1dt,1);
    t2DT=mean(T2dt,1);
    t3DT=mean(T3dt,1);
    t4DT=mean(T4dt,1);
    
%     [~,pDT,~,tDT] = ttest(DT(1:nsubjs),DT(nsubjs+1:nsubjs*2));
    
   
    hcAR=mean(HCar,1);
    t1AR=mean(T1ar,1);
    t2AR=mean(T2ar,1);
    t3AR=mean(T3ar,1);
    t4AR=mean(T4ar,1);
%     [~,pAR,~,tAR] = ttest(AR(1:nsubjs),AR(nsubjs+1:nsubjs*2));
 %% Save
    
    clusters=char(clusterNames);
    save(fullfile(savedir,['Violin5Data_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'T1fo', 'T2fo', 'T3fo', 'T4fo', 'HCfo', 'T1dt', 'T2dt', 'T3dt', 'T4dt', 'HCdt', 'T1ar', 'T2ar', 'T3ar', 'T4ar', 'HCar', 'clusters','hcDT','hcAR', 't1DT', 't2DT', 't3DT', 't4DT', 't1AR', 't2AR', 't3AR', 't4AR');%,'pDT','pAR','tDT','tAR')

end