%This script contains code to pupoulate the matrices FP/FB automatically

%need parents matrix and Kx, kr, S values
    %pre-find the combined quantities
    %find junctions, id type, set up equations

    
parents=[0 1 1 2 2]';

nL=size(parents,1);

kx=repmat(2,[nL 1]);
kr=repmat(1,[nL 1]);
L=repmat(1,[nL 1]);

a=sqrt(kr./kx);

Kr=kr.*L;
Kx=kx./L;
K=cat(1,Kx,Kr);


%% Ohm's law
iL=(1:nL)';
nC=parents>0;

vc=-cat(1,parents==0,zeros(nL,1));
dK=diag(K);
IMO=diag(ones(2*nL,1));
IMO(iL,nL+iL)=IMO(iL,nL+iL)-diag(ones(nL,1));
IMO(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;
IMCO=cat(2,vc,IMO');
C=IMO*dK*IMCO;

C1=C(iL,1);
C2=C(iL,1+iL);
C3=C(iL,iL+nL+1);

C4=diag(Kr)*(eye(5)+inv(C2)*C3);
C5=diag(Kr)*inv(C2)*C1;

%make this into fn: takes parents, hydro, gives C1,C2,C3 ... alternatively,
%C4, C5, 
%C6, Kcomp, C7...
%maybe upscaling function separate...


%% Porous pipe, point-sol'n

iL=(1:nL)';
nC=parents>0;

vc=-cat(1,parents==0,zeros(nL,1));

IMP=diag(ones(2*nL,1));
IMP(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;
IMP(sub2ind(2*[nL nL],parents(nC),nL+iL(nC)))=-1;

c0r=-a.*Kx.*tanh(a.*L/2);
c0x=a.*Kx./tanh(a.*L);
c1r=-a.*Kx.*tanh(a.*L/2);
c1x=-a.*Kx./sinh(a.*L);
cSr=2*a.*Kx.*tanh(a.*L/2);
cSx=-a.*Kx.*tanh(a.*L/2);

sFP=[2*nL 2*nL+1];

FP=zeros(sFP);
FP(sub2ind(sFP,iL,parents+1))=c1x;
FP(sub2ind(sFP,iL,iL+1))=c0x;
FP(sub2ind(sFP,iL,iL+nL+1))=cSx;
FP(sub2ind(sFP,nL+iL,parents+1))=c1r;
FP(sub2ind(sFP,nL+iL,iL+1))=c0r;
FP(sub2ind(sFP,nL+iL,iL+nL+1))=cSr;

CP=IMP*FP;

C1=C(iL,1);
C2=C(iL,1+iL);
C3=C(iL,iL+nL+1);

CL1=CP(6:10,1);
CL2=CP(6:10,2:6); %% == CP3'
CL3=CP(6:10,7:11); %% == diag(cSr)

C4=CL2*-inv(C2)*C3+CL3;
C5=CL2*-inv(C2)*C1+CL1;


%% Porous pipe, mean sol'n











