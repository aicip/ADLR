clear

nPtsPerClust = 300;
nClust = 4; % define 4 clusters
totalNumPts = nPtsPerClust*nClust;
m(:,1) = [0.1,0.2]';
m(:,2) = [0.2,0.6]';
m(:,3) = [0.8,0.9]';
m(:,4) = [0.5,0.1]';
var = .1;
bandwidth = .2;
clustMed = [];

nPtsPerClust = [100,125,150,200];
%rng(0); % comment out if you want to make the data random
x = var*randn(2,sum(nPtsPerClust));
%*** build the point set
ox = zeros(2,0);
for i = 1:nClust
    ox = [ox,repmat(m(:,i),1,nPtsPerClust(i))];   
end
x = ox + x;

% start clustering
tic
[clustCent1,point2cluster1,clustMembsCell1] = HGMeanShiftCluster(x,bandwidth,'gaussian');
toc

dict = [];
for idx = 1:size(clustCent1,2)
    a0 = clustCent1(:,idx);
%     dict = [dict a0];
    if size(clustMembsCell1{idx},2) ~=0
        merr = [];
        for i=1:size(clustMembsCell1{idx},2)
            t = clustMembsCell1{idx}(i);
            tmp = (x(:,t)-a0)'*(x(:,t)-a0);
            merr = [merr tmp];
        end
        [tmp idxsort] = sort(merr);
        a2 = x(:,clustMembsCell1{idx}(idxsort(end)));
    end
    dict = [dict a2];
end

figure;
numClust1 = length(clustMembsCell1);
hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';
% for k = 1:min(numClust1,length(cVec))
%     myMembers1 = clustMembsCell1{k};
%     myClustCen1 = clustCent1(:,k);
%     plot(x(1,myMembers1),x(2,myMembers1),[cVec(k) '>'],'MarkerSize',3)
%     plot(myClustCen1(1),myClustCen1(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
%     plot(dict(1,k),dict(2,k),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
%     legend('cluster 1','diction')
% end

for k = 1:min(numClust1,length(cVec))
    myMembers1 = clustMembsCell1{k};
    myClustCen1 = clustCent1(:,k);
    plot(x(1,myMembers1),x(2,myMembers1),[cVec(k) '>'],'MarkerSize',3)
    plot(myClustCen1(1),myClustCen1(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
    plot(dict(1,k),dict(2,k),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
    legend('cluster 1','dictionary atoms 1','dictionary atoms 1',...
            'cluster 2','dictionary atoms 2','dictionary atoms 2',...
            'cluster 3','dictionary atoms 3','dictionary atoms 3',...
            'cluster 4','dictionary atoms 4','dictionary atoms 4')
end

set(gca,'fontsize',20)




% title(['no shifting, numClust:' int2str(numClust1)])
hold off
