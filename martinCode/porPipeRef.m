function [Hx,Qx,Qr]=porPipeRef(parents,kr,kx,L,nL,Hs,Hc)

    M=findMag(parents,nL);
    
    arch=[parents L M];
    hydro=[kr kx Hs];
    a=sqrt(kr./kx);
    
    CC=linSysLF2(arch,hydro,Hc);

    Hx=Hs+CC(:,1)+CC(:,2);
    Qx=-kx.*a.*(CC(:,1)-CC(:,2));
    
    H=cat(1,Hc,Hx);
    
    c5=tanh(a.*L/2)./(a.*L);
    psiBar=c5.*(H(parents+1)+Hx)+(1-2*c5).*Hs;

    Qr=kr.*L.*(Hs-psiBar);
end