function [sinkout deltasinkout] = sinkterm(K,IM,IMT,nrootseg,psis,psicolar)

C=IM*K*IMT;
C1=C(1:nrootseg,1);
C2=C(1:nrootseg,2:nrootseg+1);
C3=C(1:nrootseg,nrootseg+2:2*nrootseg+1);
C4=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(eye(nrootseg)+...
inv(C2)*C3);

C5=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*inv(C2)*C1;
Krs=sum(sum(C4));
SUF=-C5/Krs;

C6=C4-Krs*SUF*SUF';

sinkout=C6*(psis-ones(nrootseg,1)*SUF'*psis)+Krs.*SUF.*(SUF'*psis-psicolar);
deltasinkout=C6*(psis-ones(nrootseg,1)*SUF'*psis);



end
