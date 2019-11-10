



parents=[0 1 2 3 4]';
nL=size(parents,1);

kx=repmat(1e-5,[nL 1]);
kr=repmat(1e-4,[nL 1]);
L=repmat(0.1,[nL 1]);

%% Set up matrices
[CO,CO1,CO2,CO3,CO4,CO5]=popJCOhm(parents,kr,kx,L,nL);

figure; spy(CO)


colX=[0 0 0.8];
colR=[0.8 0 0];

colPX=[0.8 0.4 0];
colPS=[0 0.65 0];

colJ=[0 100 240]/255;
colC=[240 0 160]/255;




close all
figure(1);  set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 10 8]);

hold on
set(gca,'XLim',[0 16],'YLim',[-11 1])
plot(1,-1,'s','Color',colX,'MarkerFaceColor',colX)
[i,j]=find(CO2);
plot(1+j,-i,'s','Color',colX,'MarkerFaceColor',colX)
plot(1+(1:5),-(1:5)+0.15,'s','Color',colR,'MarkerFaceColor',colR)
plot(1+(1:5),-(1:5)-0.15,'s','Color',colX,'MarkerFaceColor',colX)
[i,j]=find(CO3);
plot(6+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
[i,j]=find(CO(6:end,:));
plot(j,-5-i,'s','Color',colR,'MarkerFaceColor',colR)


line([0.5 0.5],[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.5 1],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.5 1],-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([11.5 11.5],[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([11 11.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([11 11.5],-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

plot(13-0.5,0,'s','Color',colC,'MarkerFaceColor',colC)
plot(13-0.5,-(1:5),'s','Color',colPX,'MarkerFaceColor',colPX)
plot(13-0.5,-5-(1:5),'s','Color',colPS,'MarkerFaceColor',colPS)

line([12.5 12.5]-0.5,[-10.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]-0.5,[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]-0.5,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([13.5 13.5]-0.5,[-10.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]-0.5,[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]-0.5,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

text(13.25,-5,'=','FontWeight','Bold','FontSize',18)

line([12.5 12.5]+2,[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]+2,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]+2,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([13.5 13.5]+2,[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]+2,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]+2,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

text(14.8,-1,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-2,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-3,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-4,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-5,'0','FontSize',12,'Interpreter','Latex')

plot(15,-5-(1:5),'s','Color',colJ,'MarkerFaceColor',colJ)

set(gca,'Visible','off');
set(1,'PaperUnits','centimeters','PaperPosition',[0 0 10 8])
print -depsc ~/restoredDocs/RSAstencil/images/simpleJan.eps


%% Now divide into C1,C2,C3

close all
figure(2);  set(2,'Color',[1 1 1],'Units','centimeters','Position',[5 15 10 8]);

hold on
set(gca,'XLim',[0 16],'YLim',[-11 1])
plot(1,-1,'s','Color',colX,'MarkerFaceColor',colX)
[i,j]=find(CO2);
plot(1+j,-i,'s','Color',colX,'MarkerFaceColor',colX)
plot(1+(1:5),-(1:5)+0.15,'s','Color',colR,'MarkerFaceColor',colR)
plot(1+(1:5),-(1:5)-0.15,'s','Color',colX,'MarkerFaceColor',colX)
[i,j]=find(CO3);
plot(6+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
[i,j]=find(CO(6:end,:));
plot(j,-5-i,'s','Color',colR,'MarkerFaceColor',colR)


line([0.5 0.5],[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.5 1],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.5 1],-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([11.5 11.5],[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([11 11.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([11 11.5],-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

plot(13-0.5,0,'s','Color',colC,'MarkerFaceColor',colC)
plot(13-0.5,-(1:5),'s','Color',colPX,'MarkerFaceColor',colPX)
plot(13-0.5,-5-(1:5),'s','Color',colPS,'MarkerFaceColor',colPS)

line([12.5 12.5]-0.5,[-10.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]-0.5,[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]-0.5,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([13.5 13.5]-0.5,[-10.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]-0.5,[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]-0.5,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

text(13.25,-5,'=','FontWeight','Bold','FontSize',18)

line([12.5 12.5]+2,[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]+2,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([12.5 12.75]+2,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

line([13.5 13.5]+2,[-10.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]+2,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([13.25 13.5]+2,-[10.5 10.5],'Color',[0 0 0],'LineWidth',1)

text(14.8,-1,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-2,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-3,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-4,'0','FontSize',12,'Interpreter','Latex')
text(14.8,-5,'0','FontSize',12,'Interpreter','Latex')

plot(15,-5-(1:5),'s','Color',colJ,'MarkerFaceColor',colJ)

line([1.6 1.6],[-5.4 -0.6],'Color',colX,'LineWidth',0.25)
line([1.6 6.4],[-5.4 -5.4],'Color',colX,'LineWidth',0.25)
line([1.6 6.4],[-0.6 -0.6],'Color',colX,'LineWidth',0.25)
line([6.4 6.4],[-5.4 -0.6],'Color',colX,'LineWidth',0.25)

line([0.6 0.6],[-5.4 -0.6],'Color',colC,'LineWidth',0.25)
line([0.6 1.4],[-5.4 -5.4],'Color',colC,'LineWidth',0.25)
line([0.6 1.4],[-0.6 -0.6],'Color',colC,'LineWidth',0.25)
line([1.4 1.4],[-5.4 -0.6],'Color',colC,'LineWidth',0.25)

line(5+[1.6 1.6],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)
line(5+[1.6 6.4],[-5.4 -5.4],'Color',colR,'LineWidth',0.25)
line(5+[1.6 6.4],[-0.6 -0.6],'Color',colR,'LineWidth',0.25)
line(5+[6.4 6.4],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)

line([1.6 1.6],[-5.4 -0.6]-5,'Color',colPS,'LineWidth',0.25)
line([1.6 6.4],[-5.4 -5.4]-5,'Color',colPS,'LineWidth',0.25)
line([1.6 6.4],[-0.6 -0.6]-5,'Color',colPS,'LineWidth',0.25)
line([6.4 6.4],[-5.4 -0.6]-5,'Color',colPS,'LineWidth',0.25)

set(gca,'Visible','off');
set(2,'PaperUnits','centimeters','PaperPosition',[0 0 10 8])
print -depsc ~/restoredDocs/RSAstencil/images/divideJan.eps

%% eqn for Hx, using C1,C2,C3

close all
figure(3);  set(3,'Color',[1 1 1],'Units','centimeters','Position',[5 15 16 8]);

hold on
set(gca,'XLim',[-1 20],'YLim',[-6 0])
plot(0,-(1:5),'s','Color',colPX,'MarkerFaceColor',colPX)
line(-[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(1,-3,'= -','FontWeight','Bold','FontSize',18)

line([1.6 1.6]+2,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+2,[-5.4 -5.4],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+2,[-0.6 -0.6],'Color',colX,'LineWidth',0.25)
line([6.4 6.4]+2,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)

[i,j]=find(CO2);
plot(3+j,-i,'s','Color',colX,'MarkerFaceColor',colX)
plot(3+(1:5),-(1:5)+0.15,'s','Color',colR,'MarkerFaceColor',colR)
plot(3+(1:5),-(1:5)-0.15,'s','Color',colX,'MarkerFaceColor',colX)

line(-[0.5 0.5]+4,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+4,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+4,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+8,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+8,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+8,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(8.75,-0.75,'-1','FontSize',12,'Interpreter','Latex')

line(-[0.5 0.5]+10,[-5.75 -0.25],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+10,-[0.25 0.25],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+10,-[5.75 5.75],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(CO3);
plot(9.25+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
line(-[0.5 0.5]+10.25,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+10.25,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+10.25,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+14.25,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+14.25,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+14.25,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(8.25+[1.6 1.6],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)
line(8.25+[1.6 6.4],[-5.4 -5.4],'Color',colR,'LineWidth',0.25)
line(8.25+[1.6 6.4],[-0.6 -0.6],'Color',colR,'LineWidth',0.25)
line(8.25+[6.4 6.4],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)


plot(15.5,-(1:5),'s','Color',colPS,'MarkerFaceColor',colPS)
line(15.5-[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(15.5-[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(15.5-[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+15.5,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+15.5,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+15.5,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(16.25,-3,'+','FontWeight','Bold','FontSize',18)

plot(17.5,-3,'s','Color',colC,'MarkerFaceColor',colC)

plot(18.5,-1,'s','Color',colX,'MarkerFaceColor',colX)
line(18.5-[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(18.5-[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(18.5-[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+18.5,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18.5,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18.5,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.6 0.6]+17.5,[-5.4 -0.6],'Color',colC,'LineWidth',0.25)
line([0.6 1.4]+17.5,[-5.4 -5.4],'Color',colC,'LineWidth',0.25)
line([0.6 1.4]+17.5,[-0.6 -0.6],'Color',colC,'LineWidth',0.25)
line([1.4 1.4]+17.5,[-5.4 -0.6],'Color',colC,'LineWidth',0.25)

line([0.5 0.5]+18.75,[-5.75 -0.25],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18.75,-[0.25 0.25],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18.75,-[5.75 5.75],'Color',[0 0 0],'LineWidth',1)

set(gca,'Visible','off');
set(3,'PaperUnits','centimeters','PaperPosition',[0 0 16 8])
print -depsc ~/restoredDocs/RSAstencil/images/sol123Jan.eps


%% re-combine to make C4,C5


close all
figure(4);  set(4,'Color',[1 1 1],'Units','centimeters','Position',[5 15 25 8]);

hold on
set(gca,'XLim',[0 33],'YLim',[-6 0])

line([0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5));
plot(j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(5+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(6,-3,'=','FontWeight','Bold','FontSize',18)

[i,j]=find(CO3);
plot(7+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
line(-[0.5 0.5]+8,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+8,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+8,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+12,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+12,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+12,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(6+[1.6 1.6],[-5.4 -0.6],'Color',colPS,'LineWidth',0.25)
line(6+[1.6 6.4],[-5.4 -5.4],'Color',colPS,'LineWidth',0.25)
line(6+[1.6 6.4],[-0.6 -0.6],'Color',colPS,'LineWidth',0.25)
line(6+[6.4 6.4],[-5.4 -0.6],'Color',colPS,'LineWidth',0.25)

line(-[0.5 0.5]+13.5,[-5.75 -0.25],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+13.5,-[0.25 0.25],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+13.5,-[5.75 5.75],'Color',[0 0 0],'LineWidth',1)


line(-[0.5 0.5]+14,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+14,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+14,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(13.75,-1,'1','FontSize',12,'Interpreter','Latex')
text(14.75,-2,'1','FontSize',12,'Interpreter','Latex')
text(15.75,-3,'1','FontSize',12,'Interpreter','Latex')
text(16.75,-4,'1','FontSize',12,'Interpreter','Latex')
text(17.75,-5,'1','FontSize',12,'Interpreter','Latex')

line([0.5 0.5]+18,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+18,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(19,-3,'+','FontWeight','Bold','FontSize',18)

line([1.6 1.6]+19,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+19,[-5.4 -5.4],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+19,[-0.6 -0.6],'Color',colX,'LineWidth',0.25)
line([6.4 6.4]+19,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)

[i,j]=find(CO2);
plot(20+j,-i,'s','Color',colX,'MarkerFaceColor',colX)
plot(20+(1:5),-(1:5)+0.15,'s','Color',colR,'MarkerFaceColor',colR)
plot(20+(1:5),-(1:5)-0.15,'s','Color',colX,'MarkerFaceColor',colX)

line(-[0.5 0.5]+21,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+21,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+21,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+25,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+25,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+25,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(25.75,-0.75,'-1','FontSize',12,'Interpreter','Latex')

[i,j]=find(CO3);
plot(26+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
line(-[0.5 0.5]+27,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+27,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+27,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+31,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+31,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+31,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(25+[1.6 1.6],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)
line(25+[1.6 6.4],[-5.4 -5.4],'Color',colR,'LineWidth',0.25)
line(25+[1.6 6.4],[-0.6 -0.6],'Color',colR,'LineWidth',0.25)
line(25+[6.4 6.4],[-5.4 -0.6],'Color',colR,'LineWidth',0.25)

line([0.5 0.5]+31.25,[-5.75 -0.25],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+31.25,-[0.25 0.25],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+31.25,-[5.75 5.75],'Color',[0 0 0],'LineWidth',1)

set(gca,'Visible','off');
set(4,'PaperUnits','centimeters','PaperPosition',[0 0 25 8])
print -depsc ~/restoredDocs/RSAstencil/images/c4Jan.eps

close all
figure(5);  set(5,'Color',[1 1 1],'Units','centimeters','Position',[5 15 16 8]);

hold on
set(gca,'XLim',[0 16],'YLim',[-6 0])

line([0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5,1));
plot(j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(1+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(1+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(1+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(2,-3,'=','FontWeight','Bold','FontSize',18)

[i,j]=find(CO3);
plot(2.5+j,-i,'s','Color',colR,'MarkerFaceColor',colR)
line(-[0.5 0.5]+3.5,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+3.5,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+3.5,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+7.5,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+7.5,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+7.5,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(1.5+[1.6 1.6],[-5.4 -0.6],'Color',colPS,'LineWidth',0.25)
line(1.5+[1.6 6.4],[-5.4 -5.4],'Color',colPS,'LineWidth',0.25)
line(1.5+[1.6 6.4],[-0.6 -0.6],'Color',colPS,'LineWidth',0.25)
line(1.5+[6.4 6.4],[-5.4 -0.6],'Color',colPS,'LineWidth',0.25)

line([1.6 1.6]+7,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+7,[-5.4 -5.4],'Color',colX,'LineWidth',0.25)
line([1.6 6.4]+7,[-0.6 -0.6],'Color',colX,'LineWidth',0.25)
line([6.4 6.4]+7,[-5.4 -0.6],'Color',colX,'LineWidth',0.25)

[i,j]=find(CO2);
plot(8+j,-i,'s','Color',colX,'MarkerFaceColor',colX)
plot(8+(1:5),-(1:5)+0.15,'s','Color',colR,'MarkerFaceColor',colR)
plot(8+(1:5),-(1:5)-0.15,'s','Color',colX,'MarkerFaceColor',colX)

line(-[0.5 0.5]+9,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+9,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(-[0.25 0.5]+9,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.5 0.5]+13,[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+13,-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.25 0.5]+13,-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(13.75,-0.75,'-1','FontSize',12,'Interpreter','Latex')


line(14+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(14+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(14+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

plot(15,-1,'s','Color',colX,'MarkerFaceColor',colX)

line(15+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(15+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(15+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line([0.6 0.6]+14,[-5.4 -0.6],'Color',colC,'LineWidth',0.25)
line([0.6 1.4]+14,[-5.4 -5.4],'Color',colC,'LineWidth',0.25)
line([0.6 1.4]+14,[-0.6 -0.6],'Color',colC,'LineWidth',0.25)
line([1.4 1.4]+14,[-5.4 -0.6],'Color',colC,'LineWidth',0.25)

set(gca,'Visible','off');
set(5,'PaperUnits','centimeters','PaperPosition',[0 0 15 8])
print -depsc ~/restoredDocs/RSAstencil/images/c5Jan.eps

%% Show flow calculation with C4,C5


close all
figure(6);  set(6,'Color',[1 1 1],'Units','centimeters','Position',[5 15 13 8]);

hold on
set(gca,'XLim',[0 15],'YLim',[-6 0])

line([0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5));
plot(j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(5+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(5.5+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(5.5+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(5.5+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5,1));
plot(5.5+j,-i,'s','Color',colPS,'MarkerFaceColor',colPS)

line(6.5+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(6.5+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(6.5+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

text(7.25,-3,'+','FontWeight','Bold','FontSize',18)

line(8+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(8+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(8+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5,1));
plot(8+j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(8+j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(9+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(9+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(9+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

plot(10,-3,'s','Color',colC,'MarkerFaceColor',colC)

text(11,-3,'=','FontWeight','Bold','FontSize',18)

line(12+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(12+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(12+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5,1));
plot(12+j,-i,'s','Color',colJ,'MarkerFaceColor',colJ)

line(13+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(13+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(13+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

set(gca,'Visible','off');
set(6,'PaperUnits','centimeters','PaperPosition',[0 0 13 8])
print -depsc ~/restoredDocs/RSAstencil/images/c45solJan.eps

%% three C4s

close all
figure(7);  set(7,'Color',[1 1 1],'Units','centimeters','Position',[5 15 20 8]);

hold on
set(gca,'XLim',[0 20],'YLim',[-6 0])

line([0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line([0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5));
plot(j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(5+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(5+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

line(7+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(7+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(7+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5));
plot(7+j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(7+j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(12+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(12+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(12+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)


line(14+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(14+[0.75 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(14+[0.75 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

[i,j]=find(ones(5));
plot(14+j,-i+0.075,'s','Color',colR,'MarkerFaceColor',colR)
plot(14+j,-i-0.075,'s','Color',colX,'MarkerFaceColor',colX)

line(19+[0.5 0.5],[-5.5 -0.5],'Color',[0 0 0],'LineWidth',1)
line(19+[0.25 0.5],-[0.5 0.5],'Color',[0 0 0],'LineWidth',1)
line(19+[0.25 0.5],-[5.5 5.5],'Color',[0 0 0],'LineWidth',1)

set(gca,'Visible','off');
set(7,'PaperUnits','centimeters','PaperPosition',[0 0 20 8])
print -depsc ~/restoredDocs/RSAstencil/images/threeC4s.eps



%% C6