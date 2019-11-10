

nL=5;
    iL=(1:nL)';
    sL=nL+iL;
parents=(0:4)';
nC=parents~=0;

    Kr=ones(nL,1);


    IM=diag(ones(2*nL,1));
    IM(iL,sL)=IM(iL,sL)-diag(ones(nL,1));
    IM(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;

    FB=[-1 1 -1 0 0 0 1 -1 0 0 0;...
        0 -1 1 -1 0 0 -1 1 -1 0 0;...
        0 0 -1 1 -1 0 0 -1 1 -1 0;...
        0 0 0 -1 1 -1 0 0 -1 1 -1;...
        0 0 0 0 -1 1 0 0 0 -1 1;...
        0 -1 0 0 0 0 1 0 0 0 0;...
        0 0 -1 0 0 0 0 1 0 0 0;...
        0 0 0 -1 0 0 0 0 1 0 0;...
        0 0 0 0 -1 0 0 0 0 1 0;...
        0 0 0 0 0 -1 0 0 0 0 1];
    
    %FB tries to simulate mean b/w 2 parallel roots...
    %should adapt for when Hc is not exact?
    
    
    C=IM*FB;
    C1=C(1:nL,1);
    C2=C(1:nL,2:nL+1);
    C3=C(1:nL,(nL+1+(1:nL)));
    %4-diag C2 ..  & C3? re-do rigorously, but looks like yes
    
    C4=diag(Kr)*(eye(nL)+inv(C2)*C3);
    C5=diag(Kr)*inv(C2)*C1;

    [C6,C7,SUF,Krs,Kcomp]=sufMod(C4,C5,nL);
    
    %see if C6 identical to result of Jan upscaling process on PJ model
    
    
    %never had a general way to populate matrices for the more complex cases...
    %would probably have to establish that for this... or else keep working
    %with several case studies?
    
    %try to get from full FB to upscaled FB?
    %Qs sum
    %Hx is weighted average
    %pass the coefficients through somehow to obtain sums?
    
    %Hx come from [0.1 0.3 0.4 0.2] sums to 1
    %Qs come from adding the units (1s)
    %2 eqns per layer, but n unknowns, no?
    
    
    
    nL=10;
    iL=(1:nL)';
    sL=nL+iL;
    parents=iL-1;
    nC=parents~=0;

    IM=diag(ones(2*nL,1));
    IM(iL,sL)=IM(iL,sL)-diag(ones(nL,1));
    IM(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;
    
    F1=zeros(nL,1); F1(1)=-1;
    F2=eye(nL)+diag(-ones(nL-1,1),1)+diag(-ones(nL-1,1),-1);
    F3=eye(nL)+diag(-ones(nL-1,1),1)+diag(-ones(nL-1,1),-1);
    F4=zeros(nL,1);
    F5=-eye(nL);
    F6=eye(nL);
    FB=cat(1,cat(2,F1,F2,F3),cat(2,F4,F5,F6));
    
    C=IM*FB;
    C1=C(1:nL,1);
    C2=C(1:nL,2:nL+1);
    C3=C(1:nL,(nL+1+(1:nL)));
    
    %results in 4-diag C2, C3
    
    Kr=ones(nL,1);

    C4=diag(Kr)*(eye(nL)+inv(C2)*C3);
    C5=diag(Kr)*inv(C2)*C1;

    [C6,C7,SUF,Krs,Kcomp]=sufMod(C4,C5,nL);

    
    
    