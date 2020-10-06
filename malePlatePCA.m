clear
close all

%% Script loads or generates Tierpsy feature table for manually joined male trajectories (1 trajectory per video file),
% and for the full mating plate. It also performs PCA analysis for the extracted features both at the male trajectory and at the plate summary levels.
% Author: @serenading. July 2020. 

%% Initialise
n_nonFeatVar_male = 29;
n_nonFeatVar = 28;
addpath('../AggScreening/')
useTierpsy256 = false;
maleFeatureTablePath = '/Users/sding/OneDrive - Imperial College London/wildMating/maleFeatureTable_20200730_124002.csv';
featureTablePath = '/Users/sding/OneDrive - Imperial College London/wildMating/featureTable_20200730_135546.csv';

%% Generate or load male features table

if exist(maleFeatureTablePath)
    maleFeatureTable = readtable(maleFeatureTablePath,'Delimiter',',','preserveVariableNames',true);
else
    %% Load features matrix, correspondong filenames, and metadata 
    tierpsyFeatureTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/features_summary_tierpsy_trajectory_manual_20200730_124002.csv','Delimiter',',','preserveVariableNames',true);
    tierpsyFileTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/filenames_summary_tierpsy_trajectory_manual_20200730_124002.csv','Delimiter',',','CommentStyle','#','preserveVariableNames',true);
    metadataTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/metadata_wildMating.csv','Delimiter',',','preserveVariableNames',true);
    
    %% Join tables 
    % join the Tierpsy tables to match filenames with file_id
    combinedTierpsyTable = outerjoin(tierpsyFileTable,tierpsyFeatureTable,'MergeKeys',true);
    % get just the filenames from the full path in the combined Tierpsy table
    [~, fileNamesTierpsy] = cellfun(@fileparts, combinedTierpsyTable.filename, 'UniformOutput', false);
    combinedTierpsyTable.filename = strrep(fileNamesTierpsy,'_featuresN','.hdf5');
    % finally, join tables to get strain names for each set of features
    maleFeatureTable = outerjoin(metadataTable,combinedTierpsyTable,'MergeKeys',true);
    
    %% Get feature values for the male from each video
    % get files that have skeletonisable manual trajectories and are marked as having male only trajectories
    fileLogInd = strcmp(maleFeatureTable.is_good,'True') & maleFeatureTable.male == 1;% & featureTable.n_skeletons>=5000;
    [filenames,uniqueFileInd] = unique(maleFeatureTable.filename(fileLogInd));
    dirnames = maleFeatureTable.dirname(fileLogInd);
    dirnames = dirnames(uniqueFileInd);
    % pre-allocate
    rowInd = NaN(1,numel(filenames));
    % go through each file
    for fileCtr = 1:numel(filenames)
        filename = filenames{fileCtr};
        filename_fullpath = ['/Volumes/diskAshurDT/behavgenom/Serena/wildMating/Results/' dirnames{fileCtr} '/', strrep(filename,'.hdf5','_featuresN.hdf5')];
        % load trajData
        trajData = h5read(filename_fullpath,'/trajectories_data');
        % get the male single worm manual indices
        maleIdx = unique(trajData.worm_index_manual(trajData.worm_label==1));
        assert(numel(maleIdx) ==1, 'More than 1 male manual index found.')
        rowIdx = find(strcmp(maleFeatureTable.filename, filename) & maleFeatureTable.worm_index == maleIdx);
        % record value row index
        rowInd(fileCtr) = rowIdx;
    end
    % trim featureTable down to those of male trajectories only
    maleFeatureTable = maleFeatureTable(rowInd,:);
    
    %% save male traj feature table
    writetable(maleFeatureTable,maleFeatureTablePath);
end

%% Generate or load plate worm features table
if exist(featureTablePath)
    featureTable = readtable(featureTablePath,'Delimiter',',','preserveVariableNames',true);
else
    %% load features matrix, correspondong filenames, and metadata
    tierpsyFeatureTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/features_summary_tierpsy_plate_20200730_135546.csv','Delimiter',',','preserveVariableNames',true);
    tierpsyFileTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/filenames_summary_tierpsy_plate_20200730_135546.csv','Delimiter',',','CommentStyle','#','preserveVariableNames',true);
    metadataTable = readtable('/Users/sding/OneDrive - Imperial College London/wildMating/metadata_wildMating.csv','Delimiter',',','preserveVariableNames',true);
    
    %% Join tables
    % join the Tierpsy tables to match filenames with file_id
    combinedTierpsyTable = outerjoin(tierpsyFileTable,tierpsyFeatureTable,'MergeKeys',true);
    % get just the filenames from the full path in the combined Tierpsy table
    [~, fileNamesTierpsy] = cellfun(@fileparts, combinedTierpsyTable.filename, 'UniformOutput', false);
    combinedTierpsyTable.filename = strrep(fileNamesTierpsy,'_featuresN','.hdf5');
    % finally, join tables to get strain names for each set of features
    featureTable = outerjoin(metadataTable,combinedTierpsyTable,'MergeKeys',true);
    
    %% save feature table
    writetable(featureTable,featureTablePath);
end

%% Apply feature filtering
if useTierpsy256
    feats2keep = table2cell(readtable('/Users/sding/Documents/MATLAB/AggScreening/strainsList/Tierpsy_256_short.csv','PreserveVariableNames',true,'ReadVariableNames',false));
    feats2keep = horzcat({'strainMale'},{'strainHermaphrodite'},feats2keep);
    maleFeatureTable = maleFeatureTable(:,feats2keep);
    featureTable = featureTable(:,feats2keep);
    n_nonFeatVar_male = 2;
    n_nonFeatVar = 2;
end

%% Analyze features with PCA

% extract featureMat and strain info
% for male features
maleFeatureMat = table2array(maleFeatureTable(:,n_nonFeatVar_male+1:end));
% for all plate features
featureMat = table2array(featureTable(:,n_nonFeatVar+1:end));
% for CB4856 plate features (so we can plot in CB4856's own PC space)
CB4856LogInd = strcmp(featureTable.strainHermaphrodite,'CB4856');
CB4856FeatureTable = featureTable(CB4856LogInd,:);
CB4856FeatureMat = table2array(CB4856FeatureTable(:,n_nonFeatVar+1:end));
% for MY23 plate features (so we can plot in MY23's own PC space)
MY23LogInd = strcmp(featureTable.strainHermaphrodite,'MY23');
MY23FeatureTable = featureTable(MY23LogInd,:);
MY23FeatureMat = table2array(MY23FeatureTable(:,n_nonFeatVar+1:end));

% pre-process feature matrix for PCA
%
[maleFeatureMat,~] = preprocessFeatMat(maleFeatureMat);
n_feats_male = size(maleFeatureMat,2);
%
[featureMat,~] = preprocessFeatMat(featureMat);
n_feats = size(featureMat,2);
%
[CB4856FeatureMat,~] = preprocessFeatMat(CB4856FeatureMat);
CB4856_n_feats = size(CB4856FeatureMat,2);
%
[MY23FeatureMat,~] = preprocessFeatMat(MY23FeatureMat);
MY23_n_feats = size(MY23FeatureMat,2);

% do pca
[pc_male, score_male, ~, ~, explained_male] = pca(maleFeatureMat);
[pc, score, ~, ~, explained] = pca(featureMat);
[pc_CB4856, score_CB4856, ~, ~, explained_CB4856] = pca(CB4856FeatureMat);
[pc_MY23, score_MY23, ~, ~, explained_MY23] = pca(MY23FeatureMat);

%% Male plots
% plot first two PCs and colour by different male strains
maleStrainFig = figure; hold on
CB4856MaleLogInd = strcmp(maleFeatureTable.strainMale,'CB4856');
MY23MaleLogInd = strcmp(maleFeatureTable.strainMale,'MY23');
plot(score_male(CB4856MaleLogInd,1),score_male(CB4856MaleLogInd,2),'r.')
plot(score_male(MY23MaleLogInd,1),score_male(MY23MaleLogInd,2),'b.')
xlabel(['PC1 (' num2str(round(explained_male(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_male(2))) ')%'])
legend({'CB4856 male','MY23 male'})
title(['PCA plot with ' num2str(n_feats_male) ' features'])

% plot first two PCs and shape by different mating set up
maleCrossFig = figure; hold on
maleSelfLogInd = strcmp(maleFeatureTable.strainMale,maleFeatureTable.strainHermaphrodite);
maleCrossLogInd = ~strcmp(maleFeatureTable.strainMale,maleFeatureTable.strainHermaphrodite);
plot(score_male(maleSelfLogInd,1),score_male(maleSelfLogInd,2),'ko')
plot(score_male(maleCrossLogInd,1),score_male(maleCrossLogInd,2),'k+')
xlabel(['PC1 (' num2str(round(explained_male(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_male(2))) ')%'])
legend({'same strain','different strains'})
title(['PCA plot with ' num2str(n_feats_male) ' features'])

% plot first two PCs, colour by different male strains and shape by different mating set up
combineFig = figure; hold on
CB4856MaleSelfLogInd = CB4856MaleLogInd & maleSelfLogInd;
MY23MaleSelfLogInd = MY23MaleLogInd & maleSelfLogInd;
CB4856MaleCrossLogInd = CB4856MaleLogInd & maleCrossLogInd;
MY23MaleCrossLogInd = MY23MaleLogInd & maleCrossLogInd;
plot(score_male(CB4856MaleSelfLogInd,1),score_male(CB4856MaleSelfLogInd,2),'ro')
plot(score_male(MY23MaleSelfLogInd,1),score_male(MY23MaleSelfLogInd,2),'bo')
plot(score_male(CB4856MaleCrossLogInd,1),score_male(CB4856MaleCrossLogInd,2),'r+')
plot(score_male(MY23MaleCrossLogInd,1),score_male(MY23MaleCrossLogInd,2),'b+')
xlabel(['PC1 (' num2str(round(explained_male(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_male(2))) ')%'])
legend({'CB4856 same strain','MY23 same strain','CB4856 different strains','MY23 different strains'})
title(['PCA plot with ' num2str(n_feats_male) ' features'])
set(gca,'FontSize',15,'LineWidth',1)

% Plot first three PCs and colour by different male strains
figure;
scatter3(score_male(CB4856MaleLogInd,1),score_male(CB4856MaleLogInd,2),score_male(CB4856MaleLogInd,3),'r.')
hold on
scatter3(score_male(MY23MaleLogInd,1),score_male(MY23MaleLogInd,2),score_male(MY23MaleLogInd,3),'b.')
xlabel(['PC1 (' num2str(round(explained_male(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_male(2))) ')%'])
zlabel(['PC3 (' num2str(round(explained_male(3))) ')%'])
legend({'CB4856 male','MY23 male'})
title(['PCA plot with ' num2str(n_feats_male) ' features'])

%% Plate plots

% plot first two PCs and colour by different strains
strainFig = figure; hold on
% CB4856LogInd = strcmp(featureTable.strainHermaphrodite,'CB4856');
% MY23LogInd = strcmp(featureTable.strainHermaphrodite,'MY23');
plot(score(CB4856LogInd,1),score(CB4856LogInd,2),'r.')
plot(score(MY23LogInd,1),score(MY23LogInd,2),'b.')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'CB4856','MY23'})
title(['PCA plot with ' num2str(n_feats) ' features'])

% plot first two PCs and shape by different mating set up
crossFig = figure; hold on
selfLogInd = strcmp(featureTable.strainMale,featureTable.strainHermaphrodite);
crossLogInd = ~strcmp(featureTable.strainMale,featureTable.strainHermaphrodite);
plot(score(selfLogInd,1),score(selfLogInd,2),'ko')
plot(score(crossLogInd,1),score(crossLogInd,2),'k+')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'same strain','different strains'})
title(['PCA plot with ' num2str(n_feats) ' features'])

% plot first two PCs, colour by different strains and shape by different mating set up
combineFig = figure; hold on
CB4856SelfLogInd = CB4856LogInd & selfLogInd;
MY23SelfLogInd = MY23LogInd & selfLogInd;
CB4856CrossLogInd = CB4856LogInd & crossLogInd;
MY23CrossLogInd = MY23LogInd & crossLogInd;
plot(score(CB4856SelfLogInd,1),score(CB4856SelfLogInd,2),'ro')
plot(score(MY23SelfLogInd,1),score(MY23SelfLogInd,2),'bo')
plot(score(CB4856CrossLogInd,1),score(CB4856CrossLogInd,2),'r+')
plot(score(MY23CrossLogInd,1),score(MY23CrossLogInd,2),'b+')
xlabel(['PC1 (' num2str(round(explained(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained(2))) ')%'])
legend({'CB4856 same strain','MY23 same strain','CB4856 different strains','MY23 different strains'})
title(['PCA plot with ' num2str(n_feats) ' features'])
set(gca,'FontSize',15,'LineWidth',1)

%% Plate plots (showing selfing vs. cross mating set up in each strain's own PC space)
% CB4856 mating with males from the same vs. different strains in CB4856 plate PC space
CB4856CombineFig = figure; hold on
CB4856OwnPCSelfLogInd = strcmp(CB4856FeatureTable.strainMale,'CB4856');
CB4856OwnPCCrossLogInd = strcmp(CB4856FeatureTable.strainMale,'MY23');
plot(score_CB4856(CB4856OwnPCSelfLogInd,1),score_CB4856(CB4856OwnPCSelfLogInd,2),'ro')
plot(score_CB4856(CB4856OwnPCCrossLogInd,1),score_CB4856(CB4856OwnPCCrossLogInd,2),'r+')
xlabel(['PC1 (' num2str(round(explained_CB4856(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_CB4856(2))) ')%'])
legend({'CB4856 same strain','CB4856 different strains'})
title(['PCA plot with ' num2str(CB4856_n_feats) ' features'])
set(gca,'FontSize',15,'LineWidth',1)

% MY23 mating with males from the same vs. different strains in MY23 plate PC space
MY23CombineFig = figure; hold on
MY23OwnPCSelfLogInd = strcmp(MY23FeatureTable.strainMale,'MY23');
MY23OwnPCCrossLogInd = strcmp(MY23FeatureTable.strainMale,'CB4856');
plot(score_MY23(MY23OwnPCSelfLogInd,1),score_MY23(MY23OwnPCSelfLogInd,2),'bo')
plot(score_MY23(MY23OwnPCCrossLogInd,1),score_MY23(MY23OwnPCCrossLogInd,2),'b+')
xlabel(['PC1 (' num2str(round(explained_MY23(1))) ')%'])
ylabel(['PC2 (' num2str(round(explained_MY23(2))) ')%'])
legend({'MY23 same strain','MY23 different strains'})
title(['PCA plot with ' num2str(MY23_n_feats) ' features'])
set(gca,'FontSize',15,'LineWidth',1)

%% Check to see what's inside the first PC of the hermaphrodite traj since they separate so very well
% see what's inside the first PC
[feat,featInd] = sort(pc(:,1)); % PC1 
disp('Top 10 features inside PC1 for hermaphrodite features are:')
featureTable.Properties.VariableNames(featInd(1:10))' % display top 10 feats inside PC1
 
%% Find the frame with a particular traj number
% clc
% filename = '/Volumes/diskAshurDT/behavgenom/Serena/wildMating/Results/20180819_wildMating/wildMating5.3_MY23_self_MY23_self_PC1_Ch2_19082018_131627_featuresN.hdf5';
% trajData = h5read(filename,'/trajectories_data');
% unique(trajData.worm_index_manual(trajData.worm_label==1))
% trajData.frame_number(trajData.worm_index_manual==1451)

% % get male single worm logical index
% maleSingleLogInd = trajData.worm_label==1;
% % get mating cluster logical index
% matingClusterLogInd = trajData.worm_label==2;

%% Tracking issues:
% Overlapping worms during mating, marked as cluster.
% Resolution of tail tip for males
% Single worm male coiling/bean shape not skeletonised.

%% Results: 
% Male trajectories are present in 44/60 movies. An additional movie has male in mating state only with no isolated trajectories. 
% Mating events are prevalent in videos with male trajectories (37/44). 
% No overt difference in trajectory level Tierpsy features (full or 256) for males from either strain.
% No overt difference in trajectory level Tierpsy features (full or 256) for males mating with hermaphrodites from its own strain vs. a different strain.
% Plate level Tierpsy feautures (full or 256) for hermaphrodites separate well between the two strains.
% No overt difference in plate level Tierpsy features (full or 256) in the prescence of a male from the hermaphrodites' own strain vs. a different strain.
% The lack of difference in plate level Tierpsy features still holds even when plotting each hermaphrodite strain in its own PC space.