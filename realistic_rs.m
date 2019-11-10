%vauclc ; clear ; 

WD = pwd;

addpath(genpath(WD));

%RootSys = [WD,'/in/RootSys']; % Maize root system of Valentin
% RootSys = [WD,'/in/RootSys_axelle']; % Lupin root system of Axelle
 RootSys = [WD,'/in/RootSys_fibrous']; % Fibrous root system of Valentin

CR=[WD,'/in/CondRoot.in'];

[tip,seg,age,] = read_rootsys(RootSys);
[tip,seg] = change_order(tip,seg);

IM_old = construct_IM(seg);

[Kr,Kx,kr,kx] = condroot2K(CR,seg,max(max(seg(:,end)),age));

nrootseg=length(Kr);

seg_age = age-seg(:,end);

IM=[[IM_old(2:end,1:nrootseg),diag(sparse(-ones(nrootseg,1)))];[sparse(zeros(nrootseg,nrootseg)),diag(sparse(ones(nrootseg,1)))]];

figure
subplot(1,3,1)
draw_root(RootSys)
subplot(1,3,2)
hold on
box on
subplot(1,3,3)
hold on
box on

C={'b','k'};

for i = 1:2
    pos = find(seg(:,6)==i);
    
    kr_temp=kr(pos);
    kx_temp=kx(pos);
    
    [s_age_sorted,ord]=sort(seg_age(pos));
    subplot(1,3,2)
    plot(s_age_sorted,log10(kr_temp(ord)),'--','color',C{i})
    
    subplot(1,3,3)
    plot(s_age_sorted,log10(kx_temp(ord)),'--','color',C{i})
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
figure
plot(SUF,SUF_old,'x')

psicolar=-1;
psis1=0.0*ones(nrootseg,1);
%psis1=-0.5*ones(nrootseg,1)+0.1*rand(nrootseg,1);

psix=-(C2)\(C3*psis1+C1*psicolar);
sinkdoussan1=K(nrootseg+1:2*nrootseg,nrootseg+1:2*nrootseg)*(psis1-psix);

Krstest=sum(sinkdoussan1);

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

figure
hist(log10(kcompcouvreur_mat),100);
hold on
line(log10([Kcomp_old Kcomp_old]),get(gca,'YLim'),'Color',[1 0 0])
xlabel('k_{comp}')