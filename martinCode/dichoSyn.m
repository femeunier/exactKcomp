

    %parents=[0 1 1 2 4 4 3 7 7 5 10 10 6 13 13 8 16 16 9 19 19 11 12 16 15 17 18 20 21]'; 
    %not organised monotonically
    %failed in both porPipeRef and in popJC funcitons... likely popJCMP()
    
    %parents=[0 1 1 2 3 4 4 5 5 6 7 8 9 10 10 11 11 12 12 13 13 14 15 16 17 18 19 20 21]';
    parents=[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7];
    
    %make properties a function of M?
    nL=size(parents,1);
    kr=exp(0:0.1:0.7);
    Kr=cat(1,reshape(repmat(kr(1:7),[3 1]),[21 1]),repmat(kr(8),[8 1]));
    

    kx=0.25*exp(0.7:-0.1:0);
    Kx=cat(1,reshape(repmat(kx(1:7),[3 1]),[21 1]),repmat(kx(8),[8 1]));

    L=repmat(0.1,[nL 1]);
    
    %7x3
    %8x1

    %% upscaled version of FB...
    
    
    
    
    
    
    
    
    %% setup matrices 
    
    [CO,CO1,CO2,CO3,CO4,CO5]=popJCOhm(parents,Kr,Kx,L,nL);
    [CP,CP1,CP2,CP3,CP4,CP5,CL1,CL2,CL3]=popJCPP(parents,Kr,Kx,L,nL);
    [CB,CB1,CB2,CB3,CB4,CB5]=popJCMP(parents,Kr,Kx,L,nL);

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
    
    
    
    figure; imagesc(CB4)
    figure; imagesc(CP4)
    figure; imagesc(CO4)
    
    figure; imagesc(CP4-CO4)
    norm(CP4-CO4)
    figure; imagesc(CP4-CB4)
    norm(CP4-CB4)
    
    
    [CO6,CO7,SUFO,KrsO,KcompO]=sufMod(CO4,CO5,nL);
    [CP6,CP7,SUFP,KrsP,KcompP]=sufMod(CP4,CP5,nL);
    [CB6,CB7,SUFB,KrsB,KcompB]=sufMod(CB4,CB5,nL);
    
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
    norm(CP8-CO8) %/norm(CP8) pretty low, actually
    figure; imagesc(CP8-CB8)
    norm(CP8-CB8)
    
    
    %% 
    
    
    Hc=-0.5e6;
    Hs=-0.2e6-rand(nL,1)*0.1e6;
    
    QrO=CO4*Hs+CO5*Hc;
    QrP=CP4*Hs+CP5*Hc;
    
    norm(QrO-QrP)/norm(QrP)
    
    %%
    addpath ~/Desktop/oldWork/code/

    [HxLFR,QxLFR,QrLFR]=porPipeRef(parents,Kr,Kx,L,nL,Hs,Hc); 
    
    size(Hc)
    
    C0