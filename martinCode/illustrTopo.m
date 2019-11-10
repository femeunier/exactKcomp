function [X,Y]=illustrTopo(M,parents,nL)

    
    %[L,~]=getLens(parents,M,nL);
    %P=max(L);
    J=unique(parents(2:end));
    %divide 2pi into M proportions
    
    thS=zeros(nL,2);
    thS(1,:)=[0 pi];
    
    
    r=zeros(nL,2);
    r(1,:)=[0 1];
    th=zeros(nL,2);
    th(1,2)=pi/2;
    for i=J'
        
        d=find(parents==i);
        
        thSP=thS(i,:);
        rP=r(i,2);
        thP=th(i,2);
        
        m=M(d);
        pTH=m/sum(m);
        
        dTH=(thSP(2)-thSP(1))*pTH;
        thD1=[thSP(1) thSP(1)+dTH(1)];
        thD2=[thSP(1)+dTH(1) thSP(2)];
        thD=cat(1,thD1,thD2);
        thS(parents==i,:)=thD;
        r(parents==i,:)=repmat([rP rP+1],[2 1]);
        TH=zeros(2,1);
        for j=1:2
            if m(j)==1
                TH(j)=mean(thD(j,:));
            else
                mm=M(parents==d(j));
                p=mm(1)/sum(mm);
                TH(j)=thD(j,1)+p*(thD(j,2)-thD(j,1));
            end
        end
        
        th(parents==i,:)=cat(2,repmat(thP,[2 1]),TH);
    end
  
    [X,Y]=pol2cart(th,r);
    
    
%     h=figure; hold on;
%     for i=1:nL
%         plot(X(i,:),Y(i,:))
%     end
%     
end