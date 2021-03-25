%% Contour Stim

RspContour = cat(1,RspContourA,RspContourB);
MeanRspCon = cat(1,squeeze(mean(RspContourA,2)),squeeze(mean(RspContourB,2)));

CV = 65:128;
CN = 17:64;
Bar = 1:16;

CVSI = ( max(MeanRspCon(:,65:128)') - max(MeanRspCon(:,1:64)') ) ./ ( max(MeanRspCon(:,65:128)') + max(MeanRspCon(:,1:64)') );
CNSI = ( max(MeanRspCon(:,17:64)') - max(MeanRspCon(:,[1:16 65:128])') ) ./ ( max(MeanRspCon(:,17:64)') + max(MeanRspCon(:,[1:16 65:128])') );
CVCNI = ( max(MeanRspCon(:,65:128)') - max(MeanRspCon(:,17:64)') ) ./ ( max(MeanRspCon(:,65:128)') + max(MeanRspCon(:,17:64)') );
CVPII = ( max(MeanRspCon(:,65:128)') - max(MeanRspCon(:,129:160)') ) ./ ( max(MeanRspCon(:,65:128)') + max(MeanRspCon(:,129:160)') );

%%  Response

color = jet(100);
B = zeros(8*60,20*60,3);
for ci = 1:535
    ci
    Scale(ci) = ceil(max(MeanRspCon(ci,:))/0.2) * 0.2;
    Scale(ci) = max(Scale(ci),0.8);
    Scale(ci) = min(Scale(ci),1.6);
    for i = 1:20
        for j = 1:8
            temp = 8*(i-1) + j;
            tIm = StimContour(:,:,temp);
            tIm = double(tIm) / 256;
            tIm = 0.5 - (0.5 - tIm)* 4;
            tIm = tIm(21:80,21:80);
            R = max(round(MeanRspCon(ci,temp)*100/Scale(ci)),1);
            R = min(R,100);
            B(60*(j-1)+1:60*j,60*(i-1)+1:60*i,1) = tIm * color(R,1) * 2;
            B(60*(j-1)+1:60*j,60*(i-1)+1:60*i,2) = tIm * color(R,2) * 2;
            B(60*(j-1)+1:60*j,60*(i-1)+1:60*i,3) = tIm * color(R,3) * 2;
        end
    end
    imshow(B);
    pause;
end


%%  Clustering

for i = 1:292
    t = CCA.PixelIdxList{i};
    x = ceil(t/512);
    y = mod(t-1,512) + 1;
    CxA(i) = mean(x);
    CyA(i) = mean(y);
end

for i = 1:243
    t = CCB.PixelIdxList{i};
    x = ceil(t/512);
    y = mod(t-1,512) + 1;
    CxB(i) = mean(x);
    CyB(i) = mean(y);
end

for i = 1:292
    for j = 1:292
        PairDistanceA(i,j) = 850/512 * sqrt( (CxA(i)-CxA(j))^2 + (CyA(i)-CyA(j))^2);
        CorA(i,j) = corr(MeanRspCon(i,1:128)',MeanRspCon(j,1:128)');
        DeltaCVSIA(i,j) = abs( CVSI(i) - CVSI(j) );
        DeltaCNSIA(i,j) = abs( CNSI(i) - CNSI(j) );
    end
end

for i = 1:243
    for j = 1:243
        PairDistanceB(i,j) = 850/512 * sqrt( (CxB(i)-CxB(j))^2 + (CyB(i)-CyB(j))^2);
        CorB(i,j) = corr(MeanRspCon(i+292,1:128)',MeanRspCon(j+292,1:128)');
        DeltaCVSIB(i,j) = abs( CVSI(i+292) - CVSI(j+292) );
        DeltaCNSIB(i,j) = abs( CNSI(i+292) - CNSI(j+292) );
    end
end

%%  Tuning correlation vs Pairwise Distance

x = [DistanceA(:)' DistanceB(:)'];
y = [CorA(:)' CorB(:)'];

figure,
set(gcf,'color',[1 1 1])
hold on;
for i = 1:10
    meanY(i) = mean(y(find(x>=100*(i-1)&x<100*i)));
    steY(i) = std(y(find(x>=100*(i-1)&x<100*i)))/sqrt(length((y(find(x>=100*(i-1)&x<100*i)))));
end

errorbar(50:100:950,meanY,steY,'k-','linewidth',2)
set(gca,'box','off','linewidth',2,'xticklabel',[],'yticklabel',[],'ylim',[0 0.5],'ytick',0:0.25:0.5,'xtick',0:200:1000,'xlim',[0 1000])

%%  permutation test

rng('shuffle')

for k = 1:100000
%     k
    relabel = randperm(length(y));
    y_perm = y(relabel);
    for i = 1:10
        mean_y_perm(k,i) = mean(y_perm(find(x>=100*(i-1)&x<100*i)));
    end
end

for i = 1:10
    t = sort(mean_y_perm(:,i),'ascend');
    critical_y_high(i) = t(100);
    t = sort(mean_y_perm(:,i),'descend');
    critical_y_low(i) = t(100);
    if meanY(i) > critical_y_high(i) || meanY(i) < critical_y_low(i)
        hold on;
        plot(100*i-50,meanY(i),'r*','markersize',16,'linewidth',1.2);
    end
end


%%  CVSI vs Pairwise Distance

x = [PairDistanceA(:)' PairDistanceB(:)'];
y = [PairDeltaCVSIA(:)' PairDeltaCVSIB(:)'];

figure,
set(gcf,'color',[1 1 1])
hold on;
for i = 1:10
    meanY(i) = mean(y(find(x>=100*(i-1)&x<100*i)));
    steY(i) = std(y(find(x>=100*(i-1)&x<100*i)))/sqrt(length((y(find(x>=100*(i-1)&x<100*i)))));
end

errorbar(50:100:950,meanY,steY,'k-','linewidth',2)
set(gca,'box','off','linewidth',2,'xticklabel',[],'yticklabel',[],'ylim',[0.2 0.5],'ytick',0:0.1:1,'xtick',0:200:1000)
fc = getframe(gcf);


%%  CNSI vs Pairwise Distance


x = [PairDistanceA(:)' PairDistanceB(:)'];
y = [PairDeltaCNSIA(:)' PairDeltaCNSIB(:)'];

figure,
set(gcf,'color',[1 1 1])
hold on;
for i = 1:10
    meanY(i) = mean(y(find(x>=100*(i-1)&x<100*i)));
    steY(i) = std(y(find(x>=100*(i-1)&x<100*i)))/sqrt(length((y(find(x>=100*(i-1)&x<100*i)))));
end

errorbar(50:100:950,meanY,steY,'k-','linewidth',2)
set(gca,'box','off','linewidth',2,'xticklabel',[],'yticklabel',[],'ylim',[0.2 0.5],'ytick',0:0.1:1,'xtick',0:200:1000)



%%  CVCNI vs Bar

clear p;
for ci = 1:535
    [x(ci),ind1] = max(MeanRspCon(ci,1:16));
    [y(ci),ind2] = max(MeanRspCon(ci,17:128));
    d(ci) = x(ci) - y(ci);
    t = squeeze( [RspContour(ci,:,ind1); RspContour(ci,:,16+ind2)] / max(MeanRspCon(ci,:)) );
    p(ci) = anova1(t',[1:2],'off');
end

figure,
set(gcf,'color',[1 1 1])
scatter(x(find(p<0.05&d<0)),y(find(p<0.05&d<0)),400,'r.');
hold on;
scatter(x(find(p<0.05&d>0)),y(find(p<0.05&d>0)),400,'b.');
hold on;
scatter(x(find(p>=0.05)),y(find(p>=0.05)),400,'k.');
hold on;
plot([0 1.8],[0 1.8],'k--','linewidth',2)
set(gca,'box','off','xlim',[0 2],'ylim',[0 2],'xtick',0:1:2,'ytick',0:1:2,'linewidth',4,'xticklabel',[],'yticklabel',[])


clear p;
for ci = 1:535
    [t1,ind1] = max(MeanRspCon(ci,17:64));
    [t2,ind2] = max(MeanRspCon(ci,65:128));
    t = squeeze( [RspContour(ci,:,16+ind1); RspContour(ci,:,64+ind2)] / max(MeanRspCon(ci,:)) );
    p(ci) = anova1(t',[1:2],'off');
    d(ci) = t1 - t2;
end

CCV = find(d<0&p<0.05);
CCN = find(d>0&p<0.05);

figure,
set(gcf,'color',[1 1 1])
scatter(CVCNI(find(d<0&p<0.05)),max(MeanRspCon(find(d<0&p<0.05),1:16)')./max(MeanRspCon(find(d<0&p<0.05),:)'),200,'r.');
hold on;
scatter(CVCNI(find(d>0&p<0.05)),max(MeanRspCon(find(d>0&p<0.05),1:16)')./max(MeanRspCon(find(d>0&p<0.05),:)'),200,'b.');
hold on;
scatter(CVCNI(find(p>=0.05)),max(MeanRspCon(find(p>=0.05),1:16)')./max(MeanRspCon(find(p>=0.05),:)'),200,'k.');
hold on;
plot([-0.2 -0.2],[0 1],'k--','linewidth',2.5);
hold on;
plot([0.2 0.2],[0 1],'k--','linewidth',2.5);
set(gca,'box','off','xlim',[-1 1],'ylim',[0 1],'xtick',-1:0.5:1,'ytick',0:0.25:1,'linewidth',2,'xticklabel',[],'yticklabel',[])


%%  K-Means Stim

for i = 1:20
    R(i,:) = max(MeanRspCon(:,8*(i-1)+1:8*i)')';
end

for i = 1:535
    R(:,i) = R(:,i) / max(R(:,i));
end

evaStim = evalclusters(R,'kmeans','CalinskiHarabasz','KList',[1:10]); 
figure,
set(gcf,'color',[1 1 1]);
plot(2:10,evaStim.CriterionValues(2:end),'k','linewidth',4,'marker','.','markersize',50);
set(gca,'xlim',[2 10],'ylim',[0 6],'box','off','linewidth',4,'xtick',[2 6 10],'ytick',0:3:6,'xticklabel',[],'yticklabel',[]);

idxStim = kmeans(R,2,'maxiter',10000,'replicates',10000);

Dist = 1 - corr(R');
Y = cmdscale(Dist);
for i = 1:20
    figure,
    set(gcf,'color',[1 1 1]);
    h1 = axes('Position',[0.1 0.1 0.8 0.8],'xlim',[-1 1],'ylim',[-1 1],'xtick',-1:0.5:1,'ytick',-1:0.5:1,'linewidth',4,'xticklabel',[],'yticklabel',[]);
    tIm = zeros(100);
    if i > 2
        tIm = StimContour(:,:,8*i-1);
    else
        tIm = StimContour(:,:,8*i-3);
    end
    tIm = 2 * (128 - tIm);
    tIm = imresize(tIm,20,'bilinear') * 2;
    tIm = tIm(501:1500,501:1500);
    tIm = 255-tIm;
    T = 255 * ones(1000,1000,3,'uint8');
    if idxStim(i) == 2
        T(:,:,2) = tIm;
        T(:,:,3) = tIm;
    else
        T(:,:,2) = tIm;
        T(:,:,1) = tIm;
    end
    h2 = axes('Position',[0.40+0.4*Y(i,1) 0.40+0.4*Y(i,2) 0.2 0.2]);
    imshow(T);
    fc = getframe(gcf);
    CData(:,:,:,i) = fc.cdata;
end
close all;
figure,
set(gcf,'color',[1 1 1]);
A = min(CData,[],4);
imshow(A)

for i = 1:size(Y,2)
    i
    Dy = dist(Y(:,1:i)');
    Stress(i) = sum(sum((Dist-Dy).^2)) / sum(sum((Dist).^2)); 
end

figure,
set(gcf,'Color',[1 1 1]);
plot(Stress,'k*-','linewidth',2);
set(gca,'Box','off','linewidth',4,'ylim',[0 0.4],'xtick',0:5:15,'ytick',[0 0.2 0.4],'xticklabel',[],'yticklabel',[],'xlim',[0 15]);


%%  K-Means Cell

evaCell = evalclusters(R','kmeans','CalinskiHarabasz','KList',[1:10]); 
figure,
set(gcf,'color',[1 1 1]);
plot(2:10,evaCell.CriterionValues(2:end),'k','linewidth',4,'marker','.','markersize',50);
set(gca,'xlim',[2 10],'ylim',[0 150],'box','off','linewidth',4,'xtick',[2 6 10],'ytick',0:75:150,'xticklabel',[],'yticklabel',[]);


idxCell = kmeans(R',2,'maxiter',10000,'replicates',10000);
bw1 = 0.0 * ones(512,512,3);
for i = 1:292
    bw = zeros(512);
    bw(CCA.PixelIdxList{i}) = 1;
    if idxCell(i) == 1
        bw1(:,:,1) = bw1(:,:,1) + bw;
        bw1(:,:,2) = bw1(:,:,2) - bw;
        bw1(:,:,3) = bw1(:,:,3) - bw;
    else
        bw1(:,:,1) = bw1(:,:,1) - bw;
        bw1(:,:,2) = bw1(:,:,2) - bw;
        bw1(:,:,3) = bw1(:,:,3) + bw;
    end
end
figure,imshow(bw1)

bw1 = 0.0 * ones(512,512,3);
for i = 1:243
    bw = zeros(512);
    bw(CCB.PixelIdxList{i}) = 1;
    if idxCell(i+292) == 1
        bw1(:,:,1) = bw1(:,:,1) + bw;
        bw1(:,:,2) = bw1(:,:,2) - bw;
        bw1(:,:,3) = bw1(:,:,3) - bw;
    else
        bw1(:,:,1) = bw1(:,:,1) - bw;
        bw1(:,:,2) = bw1(:,:,2) - bw;
        bw1(:,:,3) = bw1(:,:,3) + bw;
    end
end
figure,imshow(bw1)


%%  Grating Stim

MeanRspGrat= cat(1,squeeze(mean(RspGratingA,2)),squeeze(mean(RspGratingB,2)));

CRI = ( max(MeanRspGrat(:,49:54)') - max(MeanRspGrat(:,55:60)') ) ./ ( max(MeanRspGrat(:,49:54)') + max(MeanRspGrat(:,55:60)') );
CRI_1cpd = ( max(MeanRspGrat(:,49:50)') - max(MeanRspGrat(:,55:56)') ) ./ ( max(MeanRspGrat(:,49:50)') + max(MeanRspGrat(:,55:56)') );
CRI_2cpd = ( max(MeanRspGrat(:,51:52)') - max(MeanRspGrat(:,57:58)') ) ./ ( max(MeanRspGrat(:,51:52)') + max(MeanRspGrat(:,57:58)') );
CRI_4cpd = ( max(MeanRspGrat(:,53:54)') - max(MeanRspGrat(:,59:60)') ) ./ ( max(MeanRspGrat(:,53:54)') + max(MeanRspGrat(:,59:60)') );

%%  Polar vs Cart

figure,
for ci = 1:292
    [x1(ci),ind1] = max(MeanRspGrat(ci,1:48));
    [y1(ci),ind2] = max(MeanRspGrat(ci,49:60));
    d1(ci) = x1(ci) - y1(ci);
    t = squeeze( [RspGratingA(ci,:,ind1); RspGratingA(ci,:,48+ind2)] / max(MeanRspGrat(ci,1:60)) );
    p1(ci) = anova1(t',[1:2],'off');
end
set(gcf,'color',[1 1 1])
scatter(x1(find(p1<0.05&d1<0)),y1(find(p1<0.05&d1<0)),400,'r.');
hold on;
scatter(x1(find(p1<0.05&d1>0)),y1(find(p1<0.05&d1>0)),400,'b.');
hold on;
scatter(x1(find(p1>=0.05)),y1(find(p1>=0.05)),400,'k.');
hold on;


for ci = 1:243
    [x2(ci),ind1] = max(MeanRspGrat(ci+292,1:48));
    [y2(ci),ind2] = max(MeanRspGrat(ci+292,49:60));
    d2(ci) = x2(ci) - y2(ci);
    t = squeeze( [RspGratingB(ci,:,ind1); RspGratingB(ci,:,48+ind2)] / max(MeanRspGrat(ci+292,1:60)) );
    p2(ci) = anova1(t',[1:2],'off');
end

set(gcf,'color',[1 1 1])
scatter(x2(find(p2<0.05&d2<0)),y2(find(p2<0.05&d2<0)),400,'r.');
hold on;
scatter(x2(find(p2<0.05&d2>0)),y2(find(p2<0.05&d2>0)),400,'b.');
hold on;
scatter(x2(find(p2>=0.05)),y2(find(p2>=0.05)),400,'k.');
hold on;
plot([0 1.8],[0 1.8],'k--','linewidth',2)
set(gca,'box','off','xlim',[0 2],'ylim',[0 2],'xtick',0:1:2,'ytick',0:1:2,'linewidth',4,'xticklabel',[],'yticklabel',[])

%%  Distance to domain boundaries

% Map = imread('CV16X.tif');
Boundary = bwboundaries(Map);
Boundary = Boundary{1};

BoundaryDistance = zeros(1,292);
for i = 1:292
    distance = zeros(1,size(Boundary,1));
    for j = 1:size(Boundary,1)
        distance(j) = sqrt((CxA(i)-Boundary(j,1))^2+(CyA(i)-Boundary(j,2))^2);
    end
    BoundaryDistance(i) = min(distance);
    if Map(round(Cx(i)),round(Cy(i))) > 0
        BoundaryDistance(i) = -BoundaryDistance(i);
    end
end

BARCN = [BAR CN];

clear p;
clear d;
for ci = 1:292
    [x(ci),ind1] = max(MeanRspCon(ci,CV));
    [y(ci),ind2] = max(MeanRspCon(ci,BARCN));
    d(ci) = x(ci) - y(ci);
    t = squeeze( [RspContour(ci,:,CV(ind1)); RspContour(ci,:,BARCN(ind2))] / max(MeanRspCon(ci,[CV BAR CN])) );
    p(ci) = anova1(t',[1:2],'off');
end

x = 850/512 * BoundaryDistance;
y = CVSI(1:292);

figure,
set(gcf,'color',[1 1 1]);
scatter(x(find(p<0.05&d<0)),y(find(p<0.05&d<0)),100,'b.');
hold on;
scatter(x(find(p<0.05&d>0)),y(find(p<0.05&d>0)),100,'r.');
hold on;
scatter(x(find(p>=0.05)),y(find(p>=0.05)),100,'.','markeredgecolor',[.0 .0 .0]);
hold on;
plot([-1000 1000],[0.2 0.2],'k--');
hold on;
plot([-1000 1000],[-0.2 -0.2],'k--');
hold on;
plot([0 0],[-1 -0.95],'k','linewidth',3)
set(gca,'xlim',[-300 700],'ylim',[-1 1],'yticklabel',[],'xticklabel',[])



clear meanY;
clear steY;
for i = 1:9
    t = y(find(x>100*(i-1)-250&x<100*(i-1)-150));
    meanY(i) = mean(t);
    steY(i) = std(t)/sqrt(length(t))
end
figure,
set(gcf,'color',[1 1 1]);
errorbar(-200:100:600,my,sy,'linewidth',3);
hold on;
plot([0 0],[-1 1],'k--')
set(gca,'box',0,'xlim',[-300 700],'ylim',[-0.4 0.3],'ytick',-0.3:0.1:0.3,'yticklabel',[],'xticklabel',[]);


