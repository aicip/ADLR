%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% anomaly detection
%% Written by Ying Qu
%% Reference
%% Ying Qu, Wei Wang, Rui Guo, Bulent Ayhan, Chiman Kwan, Steven Vance, and Hairong Qi. 
%% "Hyperspectral Anomaly Detection through Spectral Unmixing and Dictionary based Low Rank Decomposition", 
%% accepted by IEEE Transactions on Geoscience and Remote Sensing (TGRS) 
%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
addpath(genpath('./meanshift'));
addpath(genpath('./lrr'));
addpath(genpath('./mvcnmf'));

load data/mixed;
load data/sest;
X = mixed;
[M,N,Bands] = size(X);
input = reshape(X,M*N, Bands);
input = input';
numComp = 6;
input = input/max(input(:));

% % spectral unmixing
% tic
% [Aest,sest] = do_nmfdecomp(input,numComp,M,N);
% toc

% clustering 
% % % start clustering
tic
bandwidth = 0.1;
[clustCent1,point2cluster1,clustMembsCell1] = HGMeanShiftCluster(sest,bandwidth,'gaussian');
toc

% dictionary construction 
dict = [];
for idx = 1:size(clustCent1,2)
    a0 = clustCent1(:,idx);
    dict = [dict a0];
    if size(clustMembsCell1{idx},2) ~=0
        merr = [];
        for i=1:size(clustMembsCell1{idx},2)
            x = clustMembsCell1{idx}(i);
            tmp = (sest(:,x)-a0)'*(sest(:,x)-a0);
            merr = [merr tmp];
        end
        [tmp idxsort] = sort(merr);
        a2 = sest(:,clustMembsCell1{idx}(idxsort(end)));
    end
    dict = [dict a2];
end

lambda = 0.15;

%% Dictionary based low-rank decomposition
tic 
for k=1:length(lambda)
    [A_hat,E_hat] = inexact_alm_lrr_l1(sest,dict,lambda(k),0,0);
    anomalylayer(k,:) = sum(E_hat,1)';
    anomalyimage(:,:,k) = reshape(anomalylayer(k,:),M,N);
    figure, imagesc(anomalyimage(:,:,k));
end

toc

