clc ; clear ; 

WD = pwd;

addpath(genpath(WD));

%RootSys = [WD,'/in/RootSys'];                  % Maize root system of Valentin
%RootSys = [WD,'/in/RootSys_simple'];            % Maize root system of Valentin simplified
%RootSys = [WD,'/in/RootSys_verysimple'];        % Maize root system of Valentin simplified + pruned
%RootSys = [WD,'/in/RootSys_axelle'];          % Lupin root system of Axelle
%RootSys = [WD,'/in/RootSys_fibrous'];          % Fibrous root system of Valentin
%RootSys = [WD,'/in/RootSys_taproot'];          % Taproot root system
%RootSys = [WD,'/in/RootSys_simple_fibrous'];   % simple fibrous root system
%RootSys = [WD,'/in/RootSys_simple_fibrousmultialfa'];   % simple fibrous root system
%RootSys = [WD,'/in/RootSys_simple_fibrousvertical'];   % simple fibrous root system
%RootSys = [WD,'/in/RootSys_singleroot'];       % single root
RootSys = [WD,'/in/RootSys_singlerootcompfibrous']; %single root
%RootSys = [WD,'/in/RootSys_Sunflower_3'];       % single root



%CR=[WD,'/in/CondRoot.in'];
%CR=[WD,'/in/CondRoot_lowerkx.in'];
%CR=[WD,'/in/CondRoot_mod.in'];
CR=[WD,'/in/CondRoot_singleroot.in'];
%CR=[WD,'/in/CondRoot_singleroot_varkrkx.in'];
%CR=[WD,'/in/CondRoot_singlerootmod2.in'];


[tip,seg,age,] = read_rootsys(RootSys);
[tip,seg] = change_order(tip,seg);

IM_old = construct_IM(seg);

[Kr,Kx,kr,kx] = condroot2K(CR,seg,max(max(seg(:,end)),age));


nrootseg=length(Kr);

seg_age = age-seg(:,end);

IM=[[IM_old(2:end,1:nrootseg),diag(sparse(-ones(nrootseg,1)))];[sparse(zeros(nrootseg,nrootseg)),diag(sparse(ones(nrootseg,1)))]];

figure (1);
draw_root(RootSys,1)

figure (11)
subplot(1,2,1)
hold on
box on
subplot(1,2,2)
hold on
box on

C={'b','k'};

for i = 1:2
    pos = find(seg(:,6)==i);
    
    kr_temp=kr(pos);
    kx_temp=kx(pos);
    
    [s_age_sorted,ord]=sort(seg_age(pos));
    subplot(1,2,1)
    plot(s_age_sorted,log10(kr_temp(ord)),'--','color',C{i})
    xlabel('age (d)')
    ylabel('log10 kr')
    subplot(1,2,2)
    plot(s_age_sorted,log10(kx_temp(ord)),'--','color',C{i})
    ylabel('log10 kx')
    xlabel('age (d)')
    legend('order #1', 'order #2')
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

% C4=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(eye(nrootseg)+...
% inv(C2)*C3);
% C5=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*inv(C2)*C1;

% t2 = pinv(full(C2))*C3;
% 
% if (nrootseg > 5000)
%     N=20;
%     m=mod(nrootseg,N);
%     if (m~=1)
%         X=[1:N:(nrootseg-m+1),nrootseg];
%     else
%         X=1:N:(nrootseg-m+1);
%     end
% 
%     t_all=(zeros(nrootseg));
% 
%     for i=1:(length(X)-1)
%        
%         disp(i/(length(X)-1))
%         if i < length(X)
%            cols = X(i):(X(i+1)-1);
%         else
%            cols = X(i):(X(i+1)); 
%         end
%         
%         t=sparse(zeros(nrootseg,min(length(cols),N)));
%          
%         temp=(C2\C3(:,cols)); 
%         r=repmat((1:nrootseg)',min(length(cols),N),1);
%         c=repmat(cols',nrootseg,1);
%         t=spsubsasgn(t,r,c-(i-1)*N,full(temp(:)));
%         t_all(:,cols)=t;
%     end
%     t = (t_all);
% else
%     t=mldivide(C2,C3);   
% end

C4=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(speye(nrootseg)+t);

s = mldivide(C2,C1);
C5=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*s;
Krs=sum(sum(C4));
SUF=-C5/Krs;

[Krs_old,Kcomp_old,SUF_old,~,Kr_old,Kx_old]=calculate_Krs(CR,RootSys);
figure (2);
plot(SUF,SUF_old,'x')
title('SUF comparison')
xlabel('SUF new')
ylabel('SUF old')

% Values for maize
psicolar=-4000;
psitop=-3000;
psibottom=-0;

zbottom=floor(min(seg(:,4)));


psis1=0.0*ones(nrootseg,1);
psis1=psibottom+(seg(:,4)-zbottom)*(psitop-psibottom)/(-zbottom);
%cte1=exp(4)/(exp(4)-1);
%cte2=1/(1-exp(4));

%psis1=(psibottom-psitop)*(cte1*exp((seg(:,4)-zbottom)*4/zbottom)+cte2)+psitop;



%psis1=-0.5*ones(nrootseg,1)+0.1*rand(nrootseg,1);

psix=-(C2)\(C3*psis1+C1*psicolar);
sinkdoussan1=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(psis1-psix);

psieff=SUF'*psis1;

Krstest=sum(sinkdoussan1)/(psieff-psicolar);

sinkdoussan2=C4*psis1+C5*psicolar;

C6=C4-Krs*SUF*SUF';

sinkdoussan3=C6*(psis1-ones(nrootseg,1)*SUF'*psis1)+Krs.*SUF.*(SUF'*psis1-psicolar);

sufcomp=zeros(nrootseg);
for i=1:nrootseg
    sufcomp(i,:)=C6(i,:)/C6(i,i);
    kcompcouvreur(i,i)=C6(i,i)/SUF(i)/(1-SUF(i));
    kcompcouvreur_mat(i,1)=kcompcouvreur(i,i);
end

C7=diag(1-SUF)*sufcomp+ones(nrootseg).*SUF';

sinkcouvreur=Krs*SUF*(SUF'*psis1-psicolar)+...
    kcompcouvreur*diag(SUF)*C7*(psis1-ones(nrootseg,1)*SUF'*psis1);

figure (4);
hist(log10(kcompcouvreur_mat),100);
hold on
%line(log10([Kcomp_old Kcomp_old]),get(gca,'YLim'),'Color',[1 0 0])
line(log10([Krs Krs]),get(gca,'YLim'),'Color',[1 0 0])
xlabel('log10K_{comp}')

figure (5);
plot(sinkcouvreur,sinkdoussan1,'+')
title('test sink')
xlabel('sink new')
ylabel('sink doussan equation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upscaling
dz=2;
ndepths=abs(floor(min(seg(:,4))/dz));
C6upscale=zeros(ndepths);
C7upscale=zeros(ndepths);

SUFupscale=zeros(ndepths,1);
kcompcouvreurupscale=zeros(ndepths);
sufcompupscale=zeros(ndepths);
psiupscale=zeros(ndepths,1);
sinkupscale2=zeros(ndepths,1);

%psiupscaletop=-1;
%psiupscalebottom=1;

for i=1:ndepths
    zup=-(i-1)*dz;
    zlow=-i*dz;
    mask1=seg(:,4)>=zlow & seg(:,4)<zup;
    SUFupscale(i)=sum(SUF(mask1),'double');
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
  kcompcouvreurupscaleplot(i,1)=kcompcouvreurupscale(i,i);
end
C7upscale=diag(1-SUFupscale)*sufcompupscale+ones(ndepths).*SUFupscale';

C7upscaleimage=kcompcouvreurupscale*C7upscale/Krs;
%C7upscaleimage(C7upscale==1)=0;

depths=(-dz/2:-dz:-(ndepths-1)*dz-dz/2)';

figure (6);
subplot(2,1,1)
imagesc(C7upscaleimage)
title('C7 matrix')
colorbar
subplot(2,1,2)
semilogy(depths,kcompcouvreurupscaleplot/Krs)
xlabel('depth')
ylabel('Kcomp/Krs')

psieffupscale=SUFupscale'*psiupscale;
sinkupscale=C6upscale*(psiupscale-ones(ndepths,1)*psieffupscale)...
    +Krs.*SUFupscale.*(psieffupscale-psicolar);
sinkupscaleparallel=Krs.*SUFupscale.*(psiupscale-ones(ndepths,1)*psieffupscale)+...
    Krs.*SUFupscale.*(psieffupscale-psicolar);

deltasinkupscale=sinkupscale-SUFupscale.*Krs.*(SUFupscale'*psiupscale-psicolar);
deltasinkupscaleparallel=sinkupscaleparallel-...
Krs.*SUFupscale.*(psieffupscale-psicolar);


figure (7);
plot(sinkupscale,sinkupscale2,'+')
title('test upscaled sink term')
xlabel('upscaled sinke term from summed matrix entries')
ylabel('summed small scale sink terms')




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

Kbigrootkxfit=zeros(2*ndepths);
Kbigrootkrfit=zeros(2*ndepths);
Kbigrootnofit=zeros(2*ndepths);

% upscaled parameters for a big root model are calculated using hsoil =
% zero and hcollar =-1

Qout=Krs;
hupkrfit=-1;
hupkxfit=-1;
hcollar=-1;
%hsoil=0

%%%%%%%%%
segdz=zeros(nrootseg,1);
for i=1:nrootseg
    zseg=seg(i,4);
    upnode=seg(i,5);
    %if (i<nrootseg) & (seg(i+1,1)==seg(i,5)) then
    %    downnode=seg(i+1,1);
    %end
    if upnode>0
        segdz(i)=zseg-seg(upnode,4);
    else
        segdz(i)=zseg;
    end
end
%%%%%%%

for i=1:ndepths
    zup=-(i-1)*dz;
    zlow=-i*dz;
    mask1=seg(:,4)>=zlow & seg(:,4)<zup;
        
    totlength=sum(seg(mask1,8),'double');
    kxw=sum(abs(segdz(mask1)).*kx(mask1),'double')/totlength;
    %kxw=sum(seg(mask1,8).*kx(mask1),'double')/totlength;
    Kbigrootnofit(i,i)=kxw*sum(abs(segdz(mask1)),'double')/(zup-zlow)^2;
    %Kbigrootnofit(i,i)=kxw*totlength/(zup-zlow)^2;
    Kbigrootkrfit(i,i)=Kbigrootnofit(i,i);
    Kbigrootnofit(ndepths+i,ndepths+i)=sum(Kr(mask1),'double');
    Kbigrootkxfit(ndepths+i,ndepths+i)=Kbigrootnofit(ndepths+i,ndepths+i);
    
    
    %upscaling using Kx and calculating Kr from SUF
    hx=Qout/Kbigrootkrfit(i,i)+hupkrfit;
    Kbigrootkrfit(ndepths+i,ndepths+i)=SUFupscale(i)*Krs*hcollar/hx;
    hupkrfit=hx;
    
    %upscaling using Kr and calculating Kx form SUF
    if i==1 
    %rax=(-SUFupscale(i)*Krs/Kbigrootkxfit(ndepths+i,ndepths+i)-hupkxfit)/Qout;
    rax=1/Kbigrootnofit(i,i);
    else
    rax=(-SUFupscale(i)*Krs/Kbigrootkxfit(ndepths+i,ndepths+i)-hupkxfit)/Qout;
    end
    Kbigrootkxfit(i,i)=1/rax;
    hupkxfit=-SUFupscale(i)*Krs/Kbigrootkxfit(ndepths+i,ndepths+i);
    
    %
    Qin=Qout+SUFupscale(i)*Krs*hcollar;
    Qout=Qin;
       
end    



sufkrsbigrootnofit=sufkrs(Kbigrootnofit,IMbigroot,IMTbigroot,ndepths);
sufkrsbigrootkxfit=sufkrs(Kbigrootkxfit,IMbigroot,IMTbigroot,ndepths);
sufkrsbigrootkrfit=sufkrs(Kbigrootkrfit,IMbigroot,IMTbigroot,ndepths);



[sinkbigrootnofit,deltasinkbigrootnofit]=sinkterm(Kbigrootnofit,IMbigroot,IMTbigroot,ndepths,psiupscale,psicolar);
[sinkbigrootkxfit,deltasinkbigrootkxfit]=sinkterm(Kbigrootkxfit,IMbigroot,IMTbigroot,ndepths,psiupscale,psicolar);
[sinkbigrootkrfit,deltasinkbigrootkrfit]=sinkterm(Kbigrootkrfit,IMbigroot,IMTbigroot,ndepths,psiupscale,psicolar);


figure (8);
%subplot(2,2,1)
yyaxis left
plot(depths,sinkupscale/dz,depths,sinkupscaleparallel/dz,'c--o', ...
    depths,sinkbigrootnofit/dz,'m:x');
xlabel('depth (cm)')
ylabel('sink term (cm²/d)')
%title('bigroot no fit')

yyaxis right
plot(depths,SUFupscale,depths,sufkrsbigrootnofit(1:ndepths),'x')
ylabel('SUF')
legend('sink exact','sink parallel root','sink big root','SUF','SUF big root')

%subplot(2,2,3)
%yyaxis left
%plot(depths,sinkupscale/dz,depths,sinkupscaleparallel/dz,'o', ...
%    depths,sinkbigrootkxfit/dz,'x');
%xlabel('depth (cm)')
%ylabel('sink term (cm²/d)')
%title('bigroot kx fit')

%yyaxis right
%plot(depths,SUFupscale,depths,sufkrsbigrootkxfit(1:ndepths),'x')
%ylabel('SUF')
%legend('sink exact','sink parallel root','sink big root','SUF','SUF big root')

%subplot(2,2,4)
%yyaxis left
%plot(depths,sinkupscale/dz,depths,sinkupscaleparallel/dz,'o', ...
%    depths,sinkbigrootkrfit/dz,'x');
%xlabel('depth (cm)')
%ylabel('sink term (cm²/d)')
%title('bigroot kr fit')

%yyaxis right
%plot(depths,SUFupscale,depths,sufkrsbigrootkrfit(1:ndepths),'x')
%ylabel('SUF')
%legend('sink exact','sink parallel root','sink big root','SUF','SUF big root')


%figure (9);
%yyaxis left
%plot(depths,deltasinkupscale,depths,deltasinkupscaleparallel,'o',...
%    depths,deltasinkbigroot,'x')
%xlabel('depth')
%ylabel('compensatory uptake')

%yyaxis right
%plot(depths,deltasinkupscale-deltasinkupscaleparallel,...
%    depths,deltasinkupscale-deltasinkbigroot)
%ylabel('difference compensatory uptake between models')
%legend('exact','parallel','big root','dif exact-parallel','dif exact-bigroot')

kxnofitplot=(1:1:ndepths)*0;
kxfitplot=(1:1:ndepths)*0;
krnofitplot=(1:1:ndepths)*0;
krfitplot=(1:1:ndepths)*0;

for i=1:ndepths
    kxnofitplot(i)=Kbigrootnofit(i,i);
    kxfitplot(i)=Kbigrootkxfit(i,i);
    krnofitplot(i)=Kbigrootnofit(i+ndepths,i+ndepths);
    krfitplot(i)=Kbigrootkrfit(i+ndepths,i+ndepths);
end

figure (10);
subplot (2,2,1)
plot(kxnofitplot,kxfitplot,'+')
xlabel('kx no fit');
ylabel('kx fit');

subplot (2,2,2)
plot(krnofitplot,krfitplot,'+')
xlabel('kr no fit');
ylabel('kr fit');

subplot(2,2,3)
plot(depths,kxnofitplot','o',depths,kxfitplot','+')
xlabel('depth');
ylabel('kx');
legend('kxnofit','kxfit')

subplot(2,2,4)
plot(depths,krnofitplot','o',depths,krfitplot','+')
xlabel('depth');
ylabel('kr');
legend('krnofit','krfit')

