%This script contains code to illustrate an Ohm's law parallel idealisation
%of the root system

addpath '~/restoredDocs/Brodersen/gaModMeas/'

pos=[0 0 1 0.2];
rot=0;
colX=[0 0 0.8];
colR=[0.8 0 0];

colPX=[0.8 0.4 0];
colPS=[0 0.65 0];

colC=[240 0 160]/255;

wd=0.5;
dx=0.1;

close all
figure(1); set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 10 20]);
hold on

set(gca,'XLim',[0 0.1], 'YLim',[-5.5 0])
drawResist([0.01 -1 0.03 1],90,colX,wd,dx)
drawResist([0.025 -1.25 0.05 0.5],0,colR,wd,dx)

drawResist([0.01 -2 0.03 1],90,colX,wd,dx)
drawResist([0.025 -2.25 0.05 0.5],0,colR,wd,dx)

drawResist([0.01 -3 0.03 1],90,colX,wd,dx)
drawResist([0.025 -3.25 0.05 0.5],0,colR,wd,dx)

drawResist([0.01 -4 0.03 1],90,colX,wd,dx)
drawResist([0.025 -4.25 0.05 0.5],0,colR,wd,dx)

drawResist([0.01 -5 0.03 1],90,colX,wd,dx)
drawResist([0.025 -5.25 0.05 0.5],0,colR,wd,dx)

plot(0.025,-[0 1 2 3 4 5],'sk','MarkerFaceColor',colPX)
plot(0.075,-[1 2 3 4 5],'sk','MarkerFaceColor',colPS)


text(0.005,-0.5,'$$K_x$$','FontSize',16,'Interpreter','Latex','Color',colX)
text(0.005,-1.5,'$$K_x$$','FontSize',16,'Interpreter','Latex','Color',colX)
text(0.005,-2.5,'$$K_x$$','FontSize',16,'Interpreter','Latex','Color',colX)
text(0.005,-3.5,'$$K_x$$','FontSize',16,'Interpreter','Latex','Color',colX)
text(0.005,-4.5,'$$K_x$$','FontSize',16,'Interpreter','Latex','Color',colX)

text(0.048,-0.8,'$$K_r$$','FontSize',16,'Interpreter','Latex','Color',colR)
text(0.048,-1.8,'$$K_r$$','FontSize',16,'Interpreter','Latex','Color',colR)
text(0.048,-2.8,'$$K_r$$','FontSize',16,'Interpreter','Latex','Color',colR)
text(0.048,-3.8,'$$K_r$$','FontSize',16,'Interpreter','Latex','Color',colR)
text(0.048,-4.8,'$$K_r$$','FontSize',16,'Interpreter','Latex','Color',colR)

text(0.014,-0,'$$\psi_c$$','FontSize',16,'Interpreter','Latex','Color',colC)
text(0.014,-1,'$$\psi_x$$','FontSize',16,'Interpreter','Latex','Color',colPX)
text(0.014,-2,'$$\psi_x$$','FontSize',16,'Interpreter','Latex','Color',colPX)
text(0.014,-3,'$$\psi_x$$','FontSize',16,'Interpreter','Latex','Color',colPX)
text(0.014,-4,'$$\psi_x$$','FontSize',16,'Interpreter','Latex','Color',colPX)
text(0.014,-5,'$$\psi_x$$','FontSize',16,'Interpreter','Latex','Color',colPX)

text(0.08,-1,'$$\psi_s$$','FontSize',16,'Interpreter','Latex','Color',colPS)
text(0.08,-2,'$$\psi_s$$','FontSize',16,'Interpreter','Latex','Color',colPS)
text(0.08,-3,'$$\psi_s$$','FontSize',16,'Interpreter','Latex','Color',colPS)
text(0.08,-4,'$$\psi_s$$','FontSize',16,'Interpreter','Latex','Color',colPS)
text(0.08,-5,'$$\psi_s$$','FontSize',16,'Interpreter','Latex','Color',colPS)


set(gca,'Visible','off');
set(1,'PaperUnits','centimeters','PaperPosition',[0 0 10 20])
print -depsc ~/restoredDocs/RSAstencil/images/ohmBigRoot.eps

