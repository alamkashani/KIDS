clc
clear
close all


%% step 1: determine the input data
% Code for 3-D data
% for each dataset one should change NP and -0.51+1.5*ix/NP
Datachainlink = load('dataset/chainlink.txt');
DataPoint     = Datachainlink;
%% step 2: quantize the space
NP        = 101;                %changing parameter
X         = zeros(NP,NP,NP);
Y         = zeros(NP,NP,NP);
Z         = zeros(NP,NP,NP);
D         = zeros(NP,NP,NP);
for ix = 1:NP
    for iy = 1:NP
        for iz = 1:NP
            X(ix,iy,iz) = -1.5+4*ix/NP;
            Y(ix,iy,iz) = -1.5+4*iy/NP;
            Z(ix,iy,iz) = -1.5+4*iz/NP;
        end
    end
end
X = X(:); Y = Y(:); Z = Z(:); D = D(:);
%% General Clustering
for ind =1:size(DataPoint,1)
    disp(ind);
    Data            = DataPoint(ind,:);
    rd              = sqrt((DataPoint(:,1)-Data(1)).^2+(DataPoint(:,2)-Data(2)).^2+(DataPoint(:,3)-Data(3)).^2);
    rd(ind)         =[];
    [minr,minr_ind] = min(rd);
    r               = sqrt((Data(1)-X).^2+(Data(2)-Y).^2+(Data(3)-Z).^2);
    D(r<minr)       = D(r<minr)+1/minr; % 1/minr is a changing parameter
end


D          = D./max(max(max(D)));
th         = 0.05; %changing parameter
Ima        = D;
%% cluster labeling
Ima(Ima<th)= 0;
Ima(Ima>th)= 1;
BW         = reshape(Ima,101,101,101);

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
    if Clusort(end)>(maxlength/6)
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
    Zmean= mean(Z(clusc));
    Clus_center(m,:) = [Xmean,Ymean,Zmean];
end
%% Join nearest little cluster centers to main cluster centers
CluszInd = find(ClusNumVec==0);
for n=1:numel(CluszInd)
    cluslc=Clus{CluszInd(n)};
    Xlmean= mean(X(cluslc));
    Ylmean= mean(Y(cluslc));
    Zlmean= mean(Z(cluslc));
    dis_mc = sqrt((Clus_center(:,1)-Xlmean).^2+(Clus_center(:,2)-Ylmean).^2+(Clus_center(:,3)-Zlmean).^2);
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
Z(IX) = [];
D(IX) = [];

scatter3(X,Y,Z,20,D,'filled')
Dnew(IX) = [];
figure;
scatter3(X,Y,Z,20,Dnew,'filled')
%view([1 0 0]);
figure;
scatter3(DataPoint(:,1),DataPoint(:,2),DataPoint(:,3),20,DataClusLabel,'filled')