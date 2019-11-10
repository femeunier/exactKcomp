function  [C,C1,C2,C3,C4,C5]=popJCOhm(parents,kr,kx,L,nL)

    Kr=kr.*L;
    Kx=kx./L;
    K=cat(1,Kx,Kr);

    iL=(1:nL)';
    nC=parents>0;

    vc=-cat(1,parents==0,zeros(nL,1));
    dKr=diag(Kr);
    
    IM=diag(ones(2*nL,1));
    IM(iL,nL+iL)=IM(iL,nL+iL)-diag(ones(nL,1));
    IM(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;
    IMC=cat(2,vc,IM');
    C=IM*diag(K)*IMC;

    C1=C(iL,1);
    C2=C(iL,1+iL);
    C3=C(iL,iL+nL+1);

    C4=dKr*(eye(nL)+inv(C2)*C3);
    C5=dKr*inv(C2)*C1;


end