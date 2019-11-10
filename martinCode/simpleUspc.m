

    parents=[0 1 2 2 3 4 5 5 6 6 7 8 9 10 11 11 12 12 13 13 14 14 15 16 17 18 19 20 21 22]';
    nL=size(parents,1);
    
    kr=repmat(1,[nL 1]);
    kx=repmat(0.1,[nL 1]);
    L=repmat(0.05,[nL 1]);
    

    
    %% upscaled version for MP
    
    parents=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7]';
    nL=size(parents,1);
    
    kr=repmat(1,[nL 1]);
    kx=repmat(0.1,[nL 1]);
    L=repmat(0.1,[nL 1]);

    
    
    
    %% setup matrices 
    
    [CO,CO1,CO2,CO3,CO4,CO5]=popJCOhm(parents,kr,kx,L,nL);
    [CP,CP1,CP2,CP3,CP4,CP5,CL1,CL2,CL3]=popJCPP(parents,kr,kx,L,nL);
    [CB,CB1,CB2,CB3,CB4,CB5]=popJCMP(parents,kr,kx,L,nL);

    close all
    figure(1); set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 30 10]);
    subplot(1,3,1)
    spy(CO)
    title('Discrete','FontSize',16,'Interpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    xlabel('')
    subplot(1,3,2)
    spy(CP)
    title('Point continuous','FontSize',16,'Interpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    xlabel('')
    subplot(1,3,3)
    spy(CB)
    title('Mean continuous','FontSize',16,'Interpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    xlabel('')
    set(1,'PaperUnits','centimeters','PaperPosition',[0 0 30 10])
    print -depsc ~/restoredDocs/RSAstencil/images/cComp.eps

    
    
    figure;
    subplot(1,3,1)
    spy(CO2)
    subplot(1,3,2)
    spy(CP2)
    subplot(1,3,3)
    spy(CB2)
    
    figure;
    subplot(1,3,1)
    spy(CO3)
    subplot(1,3,2)
    spy(CP3)
    subplot(1,3,3)
    spy(CB3)
    
    
    
    figure; 
    subplot(1,3,1)
    imagesc(CB4)
    subplot(1,3,2);
    imagesc(CP4)
    subplot(1,3,3);
    imagesc(CO4)
    
    close all
    figure(2); set(2,'Color',[1 1 1],'Units','centimeters','Position',[5 15 20 10]);
    subplot(1,2,1)
    imagesc(CP4-CO4)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_4$ Point-Discrete','FontSize',16,'Interpreter','Latex')
    subplot(1,2,2);
    imagesc(CP4-CB4)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_4$ Point-Mean','FontSize',16,'Interpreter','Latex')
    set(2,'PaperUnits','centimeters','PaperPosition',[0 0 20 10])
    print -depsc ~/restoredDocs/RSAstencil/images/c4Comp.eps
    
    norm(CP4-CO4)
    norm(CP4-CB4)
    
    
    [CO6,CO7,SUFO,krsO,KcompO]=sufMod(CO4,CO5,nL);
    [CP6,CP7,SUFP,krsP,KcompP]=sufMod(CP4,CP5,nL);
    [CB6,CB7,SUFB,krsB,KcompB]=sufMod(CB4,CB5,nL);
    
    figure; imagesc(CB6)
    figure; imagesc(CP6)
    figure; imagesc(CO6)
    
    figure; imagesc(CP6-CO6)
    norm(CP6-CO6)
    figure; imagesc(CP6-CB6)
    norm(CP6-CB6)
    
    figure; imagesc(CP7-CO7)
    norm(CP7-CO7)
    figure; imagesc(CP7-CB7)
    norm(CP7-CB7)
    
    CO8=diag(KcompO).*diag(SUFO)*CO7;
    CP8=diag(KcompP).*diag(SUFP)*CP7;
    CB8=diag(KcompB).*diag(SUFP)*CB7;
    
    figure; imagesc(CP8-CO8)
    norm(CP8-CO8)/norm(CP8) %pretty low, actually
    figure; imagesc(CP8-CB8)
    norm(CP8-CB8)
    
    %%
    
    close all
    figure(2); set(2,'Color',[1 1 1],'Units','centimeters','Position',[5 15 20 10]);
    subplot(1,2,1)
    imagesc(CP6-CO6)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Discrete','FontSize',16,'Interpreter','Latex')
    subplot(1,2,2);
    imagesc(CP6-CB6)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Mean','FontSize',16,'Interpreter','Latex')
    set(2,'PaperUnits','centimeters','PaperPosition',[0 0 20 10])
    print -depsc ~/restoredDocs/RSAstencil/images/c6Comp.eps
    
    norm(CP6-CO6)
    norm(CP6-CB6)

    norm(SUFP-SUFO)
    norm(SUFP-SUFB)
    
    
    
    ij={[1 2],[3 5],[4 6],[7 11],[8 12],[9 13],[10 14],[15 23],[16 24],[17 25],[18 26],[19 27],[20 28],[21 29],[22 30]}';
    [C6UP,SUFUP]=upscaleVC(ij,CP6,SUFP);
    [C6UO,SUFUO]=upscaleVC(ij,CO6,SUFO);
    [C6UB,SUFUB]=upscaleVC(ij,CB6,SUFB);
    
    close all
    figure(2); set(2,'Color',[1 1 1],'Units','centimeters','Position',[5 15 20 10]);
    subplot(1,2,1)
    imagesc(C6UP-C6UO)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Discrete','FontSize',16,'Interpreter','Latex')
    subplot(1,2,2);
    imagesc(C6UP-C6UB)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Mean','FontSize',16,'Interpreter','Latex')
    set(2,'PaperUnits','centimeters','PaperPosition',[0 0 20 10])
    print -depsc ~/restoredDocs/RSAstencil/images/c6UComp.eps
    
    norm(C6UP-C6UO)
    norm(C6UP-C6UB)

    %% 
    
    
    %% upscaled version for MP
    
    parents=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7]';
    nL=size(parents,1);
    
    kr=repmat(1,[nL 1]);
    kx=repmat(0.1,[nL 1]);
    L=repmat(0.1,[nL 1]);

    
    
    
    [CB,CB1,CB2,CB3,CB4,CB5]=popJCMP(parents,kr,kx,L,nL);

    close all
    figure(1); set(1,'Color',[1 1 1],'Units','centimeters','Position',[5 15 30 10]);
    subplot(1,3,1)
    set(gca,'Visible','off')
    subplot(1,3,2)
    set(gca,'Visible','off')
    subplot(1,3,3)
    spy(CB)
    title('Mean continuous','FontSize',16,'Interpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    xlabel('')
    set(1,'PaperUnits','centimeters','PaperPosition',[0 0 30 10])
    print -depsc ~/restoredDocs/RSAstencil/images/cUMean.eps
    
    
    [CU6,CU7,SUFU,krsU,KcompU]=sufMod(CB4,CB5,nL);
 
    close all
    figure(2); set(2,'Color',[1 1 1],'Units','centimeters','Position',[5 15 30 10]);
    subplot(1,3,1)
    imagesc(C6UP-C6UB)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Mean A','FontSize',16,'Interpreter','Latex')
    subplot(1,3,2);
    imagesc(C6UP-CU6)
    colorbar('FontSize',16,'TickLabelInterpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    title('$C_6$ Point-Mean B','FontSize',16,'Interpreter','Latex')
    subplot(1,3,3);
    spy(CB)
    title('Mean continuous','FontSize',16,'Interpreter','Latex')
    set(gca,'FontSize',16,'TickLabelInterpreter','Latex')
    xlabel('')
    set(2,'PaperUnits','centimeters','PaperPosition',[0 0 30 10])
    print -depsc ~/restoredDocs/RSAstencil/images/c6UUComp.eps

    
    
    %%
    Hc=-0.5e6;
    Hs=-0.2e6-rand(nL,1)*0.1e6;
    
    QrO=CO4*Hs+CO5*Hc;
    QrP=CP4*Hs+CP5*Hc;
    
    norm(QrO-QrP)/norm(QrP)
    
    %%
    addpath ~/Desktop/oldWork/code/

    [HxLFR,QxLFR,QrLFR]=porPipeRef(parents,kr,kx,L,nL,Hs,Hc); 
    
    size(Hc)
    
    C0
    
    
    
    
    
    