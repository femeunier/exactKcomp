

%illustrate depths
addpath '~/Desktop/oldWork/code'

%default depths from SurfEx:
layDepth=-[0 0.01,0.04,0.1,0.2,0.4,0.6,0.8,1.0,1.5,2.0,2.5,3.0]';
nLines=size(layDepth,1);

nJ=63;
[P,M,parents]=formDichoTopo(nJ);

[X,Y]=illustrTopo(M,parents,nL);
nP=size(X,1);

close all
figure(1); set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 10 20]);
hold on
for i=1:nP
    plot(0.5+sign(X(i,:)).*((-X(i,:)).^2)/100,-exp(Y(i,:)/2)/11+0.1,'Color',[0 0.6 0])
end

for i=1:nLines
    line([0 1],[layDepth(i) layDepth(i)],'LineStyle','--','Color',[0 0 0])
end
set(gca,'XColor',[1 1 1],'TickLabelInterpreter','Latex','FontSize',14,...
    'YLim',[-3 0])
ylabel('Depth (m)','Interpreter','Latex','FontSize',14)

set(1,'PaperUnits','centimeters','PaperPosition',[0 0 10 20])
print -depsc ~/restoredDocs/RSAstencil/images/surfExDepths.eps




%% second figure, for illustrating what we will upscale:


layDepth=-(0:4)';
nLines=size(layDepth,1);

nJ=7;
[P,M,parents]=formDichoTopo(nJ);
nL=2*nJ+1;

[X,Y]=illustrTopo(M,parents,nL);
nP=size(X,1);

[th,r]=cart2pol(X,Y);
figure; plot(th,r)

close all
figure(1); set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 10 20]);
hold on
for i=1:nP
    plot(X(i,:),-r(i,:),'Color',[0 0.6 0])
end

for i=1:nLines
    line([-4 4],[layDepth(i) layDepth(i)],'LineStyle','--','Color',[0 0 0])
end
set(gca,'XColor',[1 1 1],'TickLabelInterpreter','Latex','FontSize',14,...
    'YLim',[-4 0],'YTick',flipud(layDepth),'YTickLabel',fliplr({'0','-20','-40','-60','-80'}))
ylabel('Depth (cm)','Interpreter','Latex','FontSize',14)

set(1,'PaperUnits','centimeters','PaperPosition',[0 0 10 20])
print -depsc ~/restoredDocs/RSAstencil/images/upscalingNetwork.eps
