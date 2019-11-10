function [DeltaSparallel,DeltaSbigroot] = calc_model_differences(RootSys,CR)

[tip,seg,age,] = read_rootsys(RootSys);
[~,seg] = change_order(tip,seg);

IM_old = construct_IM(seg);

[Kr,Kx,kr,kx] = condroot2K(CR,seg,max(max(seg(:,end)),age));

nrootseg=length(Kr);

seg_age = age-seg(:,end);

IM=[[IM_old(2:end,1:nrootseg),diag(sparse(-ones(nrootseg,1)))];[sparse(zeros(nrootseg,nrootseg)),diag(sparse(ones(nrootseg,1)))]];


C={'b','k'};

for i = 1:2
    pos = find(seg(:,6)==i);
    
    kr_temp=kr(pos);
    kx_temp=kx(pos);
    
    [s_age_sorted,ord]=sort(seg_age(pos));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Collar_connect = zeros(1,nrootseg);
Collar_connect(seg(:,5)==0)=-1;

IMT=cat(2,[Collar_connect,zeros(1,nrootseg)]',transpose(IM));

K=diag(sparse([Kx;Kr]));
C=IM*K*IMT;
C1=C(1:nrootseg,1);
C2=C(1:nrootseg,2:nrootseg+1);
C3=C(1:nrootseg,nrootseg+2:2*nrootseg+1);
t = mldivide(C2,C3);

C4=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(speye(nrootseg)+t);

s = mldivide(C2,C1);
C5=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*s;
Krs=sum(sum(C4));
SUF=-C5/Krs;

% [Krs_old,Kcomp_old,SUF_old,~,Kr_old,Kx_old]=calculate_Krs(CR,RootSys);

psicolar=-1;
psitop=-1;
psibottom=1;

zbottom=floor(min(seg(:,4)));


psis1=psibottom+(seg(:,4)-zbottom)*(psitop-psibottom)/(-zbottom);


%psis1=-0.5*ones(nrootseg,1)+0.1*rand(nrootseg,1);

% psix=-(C2)\(C3*psis1+C1*psicolar);
% sinkdoussan1=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(psis1-psix);
% 
% psieff=SUF'*psis1;


C6=C4-Krs*(SUF*SUF');

sinkdoussan3=C6*(psis1-ones(nrootseg,1)*SUF'*psis1)+Krs.*SUF.*(SUF'*psis1-psicolar);

% sufcomp=zeros(nrootseg);
% for i=1:nrootseg
%     sufcomp(i,:)=C6(i,:)/C6(i,i);
%     kcompcouvreur(i,i)=C6(i,i)/SUF(i)/(1-SUF(i));
%     kcompcouvreur_mat(i,1)=kcompcouvreur(i,i);
% end

% C7=diag(1-SUF)*sufcomp+ones(nrootseg).*SUF';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upscaling
dz=1;
ndepths=abs(floor(min(seg(:,4))/dz));
C6upscale=zeros(ndepths);
% C7upscale=zeros(ndepths);


SUFupscale=zeros(ndepths,1);
kcompcouvreurupscale=zeros(ndepths);
sufcompupscale=zeros(ndepths);
psiupscale=zeros(ndepths,1);
sinkupscale2=zeros(ndepths,1);

for i=1:ndepths
    zup=-(i-1)*dz;
    zlow=-i*dz;
    mask1=seg(:,4)>=zlow & seg(:,4)<zup;
    SUFupscale(i)=sum(SUF(mask1));
    psiupscale(i)=SUF(mask1)'*psis1(mask1)/sum(SUF(mask1));
    sinkupscale2(i)=sum(sinkdoussan3(mask1));
  for j=1:ndepths
    zup2=-(j-1)*dz;
    zlow2=-j*dz;
    mask2=seg(:,4)>=zlow2 & seg(:,4)<zup2;
    C6upscale(i,j)=sum(sum(C6(mask1,mask2)));
  end
  sufcompupscale(i,:)=C6upscale(i,:)/C6upscale(i,i);
  kcompcouvreurupscale(i,i)=C6upscale(i,i)/SUFupscale(i)/(1-SUFupscale(i));  
end
% C7upscale=diag(1-SUFupscale)*sufcompupscale+ones(ndepths).*SUFupscale';

% C7upscaleimage=C7upscale;
% C7upscaleimage(C7upscale==1)=0;

psieffupscale=SUFupscale'*psiupscale;
sinkupscale=C6upscale*(psiupscale-ones(ndepths,1)*psieffupscale)...
    +Krs.*SUFupscale.*(psieffupscale-psicolar);
sinkupscaleparallel=Krs.*SUFupscale.*(psiupscale-ones(ndepths,1)*psieffupscale)+...
    Krs.*SUFupscale.*(psieffupscale-psicolar);

deltasinkupscale=sinkupscale-SUFupscale.*Krs.*(SUFupscale'*psiupscale-psicolar);
deltasinkupscaleparallel=sinkupscaleparallel-...
Krs.*SUFupscale.*(psieffupscale-psicolar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Big root approximation for upscaled root water uptake
IMbigroot=zeros(2*ndepths);

for i=1:ndepths-1
    IMbigroot(i,i)=1;
    IMbigroot(i+ndepths,i+ndepths)=1;
    IMbigroot(i,i+1)=-1;
    IMbigroot(i,ndepths+i)=-1;
end
IMbigroot(ndepths,ndepths)=1;
IMbigroot(ndepths,2*ndepths)=-1;
IMbigroot(2*ndepths,2*ndepths)=1;

dummy=[-1;zeros(2*ndepths-1,1)];
IMTbigroot=cat(2,dummy,transpose(IMbigroot));

Kbigroot=zeros(2*ndepths);
Qout=Krs;
hup=-1;
for i=1:ndepths
    zup=-(i-1)*dz;
    zlow=-i*dz;
    mask1=seg(:,4)>=zlow & seg(:,4)<zup;
    Kbigroot(ndepths+i,ndepths+i)=sum(Kr(mask1));
    rax=(-SUFupscale(i)*Krs/Kbigroot(ndepths+i,ndepths+i)-hup)/Qout;
    Kbigroot(i,i)=1/rax;
    hup=-SUFupscale(i)*Krs/Kbigroot(ndepths+i,ndepths+i);
    Qin=Qout-SUFupscale(i)*Krs;
    Qout=Qin;
end    

[~,deltasinkbigroot]=sinkterm(Kbigroot,IMbigroot,IMTbigroot,ndepths,psiupscale,psicolar);

DeltaSparallel= deltasinkupscale-deltasinkupscaleparallel;
DeltaSbigroot = deltasinkupscale-deltasinkbigroot;


end





