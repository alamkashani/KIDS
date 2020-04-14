clc
clear all
close all


%% step 1: determine the input data
% Code for 2-D data
% for each dataset one should change NP and 0+20*ix/N
DataAggregation = load('dataset/Aggregation.txt');
DataPoint       = DataAggregation;
%% step 2: quantize the space
NP        = 400;                %changing parameter
X         = zeros(NP,NP);
Y         = zeros(NP,NP);
D         = zeros(NP,NP);
for ix = 1:NP
    for iy = 1:NP
            X(ix,iy) = 0+40*ix/NP;
            Y(ix,iy) = 0+40*iy/NP;
    end
end
X = X(:); Y = Y(:); D = D(:);
%% General Clustering (IDS)
for ind =1:size(DataPoint,1)
    disp(ind);
    Data            = DataPoint(ind,:);
    rd              = sqrt((DataPoint(:,1)-Data(1)).^2+(DataPoint(:,2)-Data(2)).^2);
    rd(ind)         = [];
    [minr,minr_ind] = min(rd);
    r               = sqrt((Data(1)-X).^2+(Data(2)-Y).^2);
    D(r<minr)       = D(r<minr)+1/minr; % 1/minr is a changing parameter
end

%% cluster labeling
D          = D./(max(max(D)));
th         = 0.1; %changing parameter
D(D<th)    = 0;
D(D>th)    = 1;
BW         = reshape(D,NP,NP);

CC         = bwconncomp(BW);
Clus  = CC.PixelIdxList;
initClusNum =CC.NumObjects;
for i=1:initClusNum
    Cluslength(i) = length(Clus{i});
end
ClusNum =1;
ClusNumVec = zeros(size(Cluslength));
[maxlength,maxInd] = max(Cluslength);
ClusNumVec(maxInd)=1;
Clusort = sort(Cluslength);
Clusort(end)=[];
while (length(Clusort))>0
    if Clusort(end)>(maxlength/10)
        ClusNum = ClusNum+1;
        ClusNumVec((Cluslength==Clusort(end)))=ClusNum;
        Clusort(end)=[];
    else
        Clusort(end)=[];
    end
end
%% finding main cluster centers
ClusInd = find(ClusNumVec>0);
for m=1:ClusNum
    clusc=Clus{ClusInd(m)};
    Xmean= mean(X(clusc));
    Ymean= mean(Y(clusc));
    Clus_center(m,:) = [Xmean,Ymean];
end
%% Join nearest little cluster centers to main cluster centers
CluszInd = find(ClusNumVec==0);
for n=1:numel(CluszInd)
    cluslc=Clus{CluszInd(n)};
    Xlmean= mean(X(cluslc));
    Ylmean= mean(Y(cluslc));
    dis_mc = sqrt((Clus_center(:,1)-Xlmean).^2+(Clus_center(:,2)-Ylmean).^2);
    [clusr,r_ind] = min(dis_mc);
    ClusNumVec(CluszInd(n))=ClusNumVec(ClusInd(r_ind));
end

Dnew = zeros(size(X));
for j = 1:numel(Clus)
    l = Clus{j};
    Dnew(l) = ClusNumVec(j);
end

%% plot
IX = (D<th)|(Dnew<th);
X(IX) = [];
Y(IX) = [];
D(IX) = [];

figure
scatter(X,Y,20,D,'filled')
Dnew(IX) = [];
figure;
scatter(X,Y,20,Dnew,'filled')
%view([1 0 0]);
figure;
plot(DataPoint(:,1),DataPoint(:,2),'o')