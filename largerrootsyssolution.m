%global psis1 psis2;

ra=0.1;
rs=1;
nrootseg=9;

% Doussan solution

IM=...
[1	-1	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0;...
0	1	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0;...
0	0	1	-1	0	0	0	0	0	0	0	-1	0	0	0	0	0	0;...
0	0	0	1	-1	0	0	0	0	0	0	0	-1	0	0	0	0	0;...
0	0	0	0	1	0	0	0	0	0	0	0	0	-1	0	0	0	0;...
0	0	0	0	0	1	-1	0	0	0	0	0	0	0	-1	0	0	0;...
0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	-1	0	0;...
0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	-1	0;...
0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	-1;...
0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0;...
0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;...
0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;...
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0;...
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0;...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0;...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0;...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0;...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1];



IMT=cat(2,[-1;0;-1;0;0;-1;0;0;0;0;0;0;0;0;0;0;0;0],transpose(IM));

K=diag([1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra...
    0.1/rs 1/rs 0.1/rs 0.1/rs 1/rs 0.1/rs 0.1/rs 0.1/rs 1/rs]);

%K=diag([1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra 1/ra...
%    1/rs 1/rs 1/rs 1/rs 1/rs 1/rs 1/rs 1/rs 1/rs]);

C=IM*K*IMT;
C1=C(1:nrootseg,1);
C2=C(1:nrootseg,2:nrootseg+1);
C3=C(1:nrootseg,nrootseg+2:2*nrootseg+1);
C4=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(eye(nrootseg)+...
inv(C2)*C3);

C5=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*inv(C2)*C1;
Krs=sum(sum(C4));
SUF=-C5/Krs;


psicolar=-1;
psis1=[-0.5;0;-0.5;0;0.5;-0.5;0;0.5;1];



psix=-inv(C2)*(C3*psis1+C1*psicolar);
sinkdoussan1=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(psis1-psix);


sinkdoussan2=C4*psis1+C5*psicolar;

C6=C4-Krs*SUF*SUF';

sinkdoussan3=C6*(psis1-ones(nrootseg,1)*SUF'*psis1)+Krs.*SUF.*(SUF'*psis1-psicolar);

sufcomp=C6*0;
for i=1:nrootseg
    sufcomp(i,:)=C6(i,:)/C6(i,i);
    kcompcouvreur(i,i)=C6(i,i)/SUF(i)/(1-SUF(i));
end
C7=diag(1-SUF)*sufcomp+ones(nrootseg).*SUF';

sinkcouvreur=Krs*SUF*(SUF'*psis1-psicolar)+...
    kcompcouvreur*diag(SUF)*C7*(psis1-ones(nrootseg,1)*SUF'*psis1);


depth1=[1,3,6];
depth2=[2,4,7];
depth3=[5,8];
depth4=[9];
tabledepth=table(depth1,depth2,depth3,depth4);

ndepths=4;
C6upscale=zeros(ndepths);
C7upscaleincorrect=zeros(ndepths);

SUFupscale=zeros(ndepths,1);
kcompcouvreurupincorrect=zeros(ndepths);
kcompcouvreurupscale=zeros(ndepths);
sufcompupscale=zeros(ndepths);
psiupscale=zeros(ndepths,1);
sinkupscale2=zeros(ndepths,1);

for i=1:ndepths
    SUFupscale(i)=sum(SUF(table2array(tabledepth(1,i))));
        psiupscale(i)=SUF(table2array(tabledepth(1,i)))'*psis1(table2array(tabledepth(1,i)))...
        /sum(SUF(table2array(tabledepth(1,i))));
    sinkupscale2(i)=sum(sinkdoussan3(table2array(tabledepth(1,i))));
  for j=1:ndepths
      C6upscale(i,j)=sum(sum(C6(table2array(tabledepth(1,i)),table2array(tabledepth(1,j)))));
      C7upscaleincorrect(i,j)=sum(sum(C7(table2array(tabledepth(1,i)),table2array(tabledepth(1,j)))));
      kcompcouvreurupincorrect(i,j)=sum(sum(kcompcouvreur(table2array(tabledepth(1,i)),...
          table2array(tabledepth(1,j)))));
  end
end

sinkupscale=C6upscale*(psiupscale-ones(ndepths,1)*SUFupscale'*psiupscale)...
    +Krs.*SUFupscale.*(SUFupscale'*psiupscale-psicolar);

deltasinkupscale=sinkupscale-SUFupscale.*Krs.*(SUFupscale'*psiupscale-psicolar);

sinkupscaleinocrrect=kcompcouvreurupincorrect*diag(SUFupscale)*C7upscaleincorrect*(psiupscale-ones(ndepths,1)*SUFupscale'*psiupscale)...
    +Krs.*SUFupscale.*(SUFupscale'*psiupscale-psicolar);


for i=1:ndepths
    sufcompupscale(i,:)=C6upscale(i,:)/C6upscale(i,i);
    kcompcouvreurupscale(i,i)=C6upscale(i,i)/SUFupscale(i)/(1-SUFupscale(i));
end
C7upscale=diag(1-SUFupscale)*sufcompupscale+ones(ndepths).*SUFupscale';

sinkcouvreurupscale=Krs*SUFupscale*(SUFupscale'*psiupscale-psicolar)+...
    kcompcouvreurupscale*diag(SUFupscale)*C7upscale*(psiupscale-ones(ndepths,1)*SUFupscale'*psiupscale);

sinkcouvreurupscaleapprox=Krs*SUFupscale*(SUFupscale'*psiupscale-psicolar)+...
    Krs*SUFupscale.*(psiupscale-SUFupscale'*psiupscale);

deltasinkcouvreurupscaleapprox=sinkcouvreurupscaleapprox-SUFupscale.*Krs.*(SUFupscale'*psiupscale-psicolar);

deltapsiupscale=psiupscale-SUFupscale'*psiupscale;



%Big root approximation for upscaled root water uptake

IMbigroot=...
[1,-1,0,0,-1,0,0,0;...
0,1,-1,0,0,-1,0,0;...
0,0,1,-1,0,0,-1,0;...
0,0,0,1,0,0,0,-1;...
0,0,0,0,1,0,0,0,;...
0,0,0,0,0,1,0,0;...
0,0,0,0,0,0,1,0;...
0,0,0,0,0,0,0,1];

IMTbigroot=cat(2,[-1;0;0;0;0;0;0;0],transpose(IMbigroot));
Kbigroot= diag([3/ra 3/ra 2/ra 1/ra...
    0.3/rs (1/rs+0.2/rs) (1/rs+0.1/rs) 1/rs]);

%Kbigroot= diag([3/ra 3/ra 2/ra 1/ra...
%    3/rs 3/rs 2/rs 1/rs]);


%sufkrsfit=@(ra1,ra2,ra3,ra4,rs1) sufkrs(Kbigroot(ra1,ra2,ra3,ra4,rs1),IMbigroot,IMTbigroot,ndepths);

%fitfun= @(parameter) (sufkrsfit(parameter(1),parameter(2),parameter(3),parameter(4),parameter(5))-...
%    cat(1,SUFupscale,Krs));

%fitparameter=lsqnonlin(fitfun,[0.1,0.1,0.1,0.1,1],[0.00001,0.00001,0.00001,0.00001,0.00001],...
%    [100,100,100,100,100]);

%Kbigrootfit=diag([1/fitparameter(1),1/fitparameter(2),1/fitparameter(3),1/fitparameter(4),...
%    1/fitparameter(5),1/fitparameter(5),1/fitparameter(5),1/fitparameter(5)]);

sufkrsbigroot=sufkrs(Kbigroot,IMbigroot,IMTbigroot,ndepths);

[sinkbigroot,deltasinkbigroot]=sinkterm(Kbigroot,IMbigroot,IMTbigroot,ndepths,psiupscale,psicolar);


figure (1);
title('kx=10, kr=0.1 (1.0 at root tips), hybrid parallel-big root')
yyaxis left;
plot(-1*(1:4),deltasinkupscale,'b',-1*(1:4),deltasinkcouvreurupscaleapprox,'c-',...
    -1*(1:4),deltasinkbigroot,'m--');
xlabel('depth');
ylabel('delta sink')
yyaxis right;
plot(-1*(1:4),psiupscale,'r');
%xticks([-1 -2 -3 -4]);
ylabel('soil water potential')
legend('exact','Couvreur approx','Big root approx')

figure (2);
title('kx=10, kr=0.1 (1.0 at root tips), hybrid parallel-big root')
yyaxis left;
plot(-1*(1:4),sinkupscale,'b',-1*(1:4),sinkcouvreurupscaleapprox,'c-',...
    -1*(1:4),sinkbigroot,'m--');
xlabel('depth');
ylabel('sink')
yyaxis right;
plot(-1*(1:4),psiupscale,'r');
%xticks([-1 -2 -3 -4]);
ylabel('soil water potential')
legend('exact','Couvreur approx','Big root approx')

figure (3);
title('kx=10, kr=0.1 (1.0 at root tips), hybrid parallel-big root')
yyaxis left;
plot(-1*(1:4),SUFupscale,'b',...
    -1*(1:4),sufkrsbigroot(1:4),'m--');
xlabel('depth');
ylabel('SUF')
yyaxis right;
plot(-1*(1:4),psiupscale,'r');
%xticks([-1 -2 -3 -4]);
ylabel('soil water potential')
legend('exact-Couvreur approx','Big root approx')