function [Krs_s,Kcomp,SUF,c,Kr,Kx,seg_info]=calculate_Krs(condroot,rootsys,Hsoil,IM,type)

% Local conductivities

K=read_condroot(condroot);

% C matrix (doussan)

[info_tip,seg_info,age]=read_rootsys(rootsys);
if nargin<5
   type='def'; 
end
if nargin <4 || IM==-1
    IM=rootsys2C(rootsys,seg_info);
end
if nargin <3 ||Hsoil==-1
    Hsoil=-150+seg_info(:,4); 
end
PHs2=Hsoil;

[info_tip,seg_info]=change_order(info_tip,seg_info,type);
if seg_info(1,end)<age/3
    seg_info(:,end)=max(0,age-seg_info(:,end));
end
nseg=size(seg_info,1);
% Kr and Kx for each segment
Kr=zeros(size(seg_info,1),1);
Kx=zeros(size(seg_info,1),1);

axes=unique(seg_info(:,6));

for i=1:length(axes)
    pos=ismember(seg_info(:,6),axes(i));
    Kr(pos)=interp1(K{i}(:,1),K{i}(:,2),seg_info(pos,11));
    Kx(pos)=interp1(K{i+3}(:,1),K{i+3}(:,2),seg_info(pos,11));
end

% Kr(seg_info(:,11)<14)=0;
Kr=Kr.*seg_info(:,9);
Kx=Kx./seg_info(:,8);


%     N=max(seg_info(:,7));
%     for i=2:N
%         pos=find(seg_info(:,7)==i);
%         Krs_lat(i,:)=[info_tip(i,8),calculate_Krs_HA((Kr(pos)),(Kx(pos)),construct_IMrecti(length(pos)))];     
%     end

im=IM(2:end,:);
c=im*(diag(-sparse([Kx;Kr])))*im';

% Hsr=-abs(100*(randn(nseg,1)-5));
Hsr=-300*ones(nseg,1);
Hsr_U=-300*ones(nseg,1);
Hcol=-1000;

% *** Calcul des matrices C et c de Doussan


% invc=inv(c);

%Krs0=-Kx(1)*Kr'*invc(:,1);
%SUF0=-Kx(1)*Kr.*invc(:,1)/Krs0;

Q_Dou=Kr.*(Hsr+c\(Kr.*Hsr+[Kx(1)*Hcol;zeros(nseg-1,1)]));
Tact_Dou=sum(Q_Dou);
Q_Dou_U=Kr.*(Hsr_U+c\(Kr.*Hsr+[Kx(1)*Hcol;zeros(nseg-1,1)]));
SUF=abs(Q_Dou_U/sum(Q_Dou_U));

Krs_s=abs(sum((Q_Dou_U)))/(Hsr_U(1)-Hcol);

%%%%%%%%%%%%%%%%%%%%%% Calcul of Kcomp --> false now
% SSF=SUF;
% SSFunstock(:,2)=SSF;
% nElm=length(SSF);
% for i=1:nElm
%     SSFunstock(i,1)=i;
%     Hs_cube2(i)=-150+seg_info(i,4);
% end
% nk=nElm;  
% SSFt=0.0001;
% SSFcum=0;
% j=0;
% 
% while(SSFcum<0.9)
%     SSFt=0.5*SSFt+0.5*sum(SSFunstock(1:nk,2))/nk;
%     SSFunstockNew=zeros(nk,2);
%     k=0;
%     for i=1:nk
%         if (SSFunstock(i,2)>=SSFt) 
%              SSFcum=SSFcum+SSFunstock(i,2)
%              j=j+1;
%              phi=

SSF=SUF;
nrec=length(SUF);
Hseq=SSF'*PHs2;
k=0;
Jcol2=zeros(1,1);
Joutr2=zeros(nrec,1);

PHr_tot2=c\(diag(-sparse(Kr))*PHs2-[Kx(1)*Hcol;zeros(nrec-1,1)]);

for irecn=1:nrec
    Joutr2(irecn) = Kr(irecn)*(PHs2(irecn)-PHr_tot2(irecn));
    Jcol2=Jcol2+Joutr2(irecn);
end

  
for i=1:nrec
    if (SSF(i)>0.5/nrec)
        PhiSUF(i)=Joutr2(i)/SSF(i)-Jcol2;
        k=k+1;
        ind(k)=i;
    else
        PhiSUF(i)=0;
    end  
end

if exist('ind')==0
    Kcomp='NaN';
else
    Kcomp=((PHs2(ind(1:k))-Hseq)'*PhiSUF(ind(1:k))')/((PHs2(ind(1:k))-Hseq)'*(PHs2(ind(1:k))-Hseq));    
    SSFcum=sum(SSF(ind(1:k)));
end

end