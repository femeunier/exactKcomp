%This script contains code to re-create the Ohm's law Jan system (Doussan solution?)
%for Valentin's simple system, from notes.

%Qs: how to upscale the LF model?


parents=[0 1 1 2 2]';

nL=size(parents,1);

kx=repmat(2,[nL 1]);
kr=repmat(1,[nL 1]);
L=repmat(1,[nL 1]);


%% boundary conditions: Hs, Hc

Hc=-4e5;
Hs=[-2e5 -1.7e5 -1.6e5 -1.5e5 -1.4e5]';

%solution for Hx?

%% reference for LF

    addpath ~/Desktop/oldWork/code/

    M=[3 2 1 1 1]';
    arch=[parents L M];
    hydro=[Kr Kx Hs];
    a=sqrt(Kr./Kx);
    
    CC=linSysLF2(arch,hydro,Hc);

    HxLF=Hs+CC(:,1)+CC(:,2);
    QxLF=Kx.*a.*(CC(:,1)-CC(:,2));
    
    psi1=HxLF(1);
    psi0=HxLF(2);
    psiS=Hs(2);
    

    c5=tanh(a(2)*L(2)/2)/(a(2)*L(2));
    psiBar2=c5*(psi1+psi0)+(1-2*c5)*psiS;
%    -Kr(2)*L(2)*(psiBar2-psiS)
    
    Qr2=a(2)*Kx(2)*tanh(a(2)*L(2)/2)*(2*psiS-psi1-psi0);
    Q02=a(2).*Kx(2).*(psi0/tanh(a(2)*L(2))-psi1/sinh(a(2)*L(2))-tanh(a(2)*L(2)/2)*psiS);

    
%% Jan's form of Ohm's law
Kr=kr.*L;
Kx=kx./L;
K=cat(1,Kx,Kr);


IMO=[1 -1 -1 0 0 -1 0 0 0 0;...
    0 1 0 -1 -1 0 -1 0 0 0;...
    0 0 1  0  0 0  0 -1 0 0;...
    0 0 0  1  0 0  0 0 -1 0;...
    0 0 0 0 1 0 0 0 0 -1; ...
    0 0 0 0 0 1 0 0 0 0;...
    0 0 0 0 0 0 1 0 0 0;...
    0 0 0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 0 0 1]; %could probably fill this from parents!

IMCO=cat(2,[-1 0 0 0 0 0 0 0 0 0]',IMO');

FO=diag(K)*IMCO;
C=IMO*FO;


C1=C(1:5,1);
C2=C(1:5,2:6);
C3=C(1:5,7:11);

diag(K)*IMCO


HxJO=-inv(C2)*(C3*Hs+Hc*C1); 

KR=Kr.*L;
C4=diag(KR)*(eye(5)+inv(C2)*C3);
C5=diag(KR)*inv(C2)*C1;
%C4*Hs+C5*Hc checks out == Qr

KrsOhm=-sum(C5);
%sum(sum(C4))


SUFohm=C5/sum(C5);
%sum(C4,2)/sum(sum(C4))
HeffOhm=SUFohm'*Hs;

C6=C4-KrsOhm*SUFohm*SUFohm';

C6*Hs+KrsOhm*(HeffOhm-Hc)*SUFohm


SUFO2=[SUFohm(1:3); sum(SUFohm(4:5))];
HeffO2=SUFohm'*Hs2;

QrO2=C6*Hs2+KrsOhm*(HeffO2-Hc)*SUFohm


KcompOhm=diag(C6)./(SUFohm.*(1-SUFohm));
C7=diag((1-SUFohm)./diag(C6))*C6+ones(5,1)*SUFohm';

[diag(KcompOhm)*diag(SUFohm)*C7*(Hs-HeffOhm)+KrsOhm*(HeffOhm-Hc)*SUFohm C4*Hs+C5*Hc]


%% and the same in porous pipe...


IMP=[1 -1 -1 0 0 0 -1 -1 0 0;...
    0 1 0 -1 -1 0 0 0 -1 -1;...
    0 0 1 0 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 1 0 0 0 0 0; ...
    0 0 0 0 0 1 0 0 0 0;...
    0 0 0 0 0 0 1 0 0 0;...
    0 0 0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 0 0 1]; %sign varies by author...

c0r=-a.*kx.*tanh(a.*L/2);
c0x=a.*kx./tanh(a.*L);
c1r=-a.*kx.*tanh(a.*L/2);
c1x=-a.*kx./sinh(a.*L);
cSr=2*a.*kx.*tanh(a.*L/2);
cSx=-a.*kx.*tanh(a.*L/2);


FP=[c1x(1) c0x(1) 0 0 0 0 cSx(1) 0 0 0 0;...
    0 c1x(2) c0x(2) 0 0 0 0 cSx(2) 0 0 0;...
    0 c1x(3) 0 c0x(3) 0 0 0 0 cSx(3) 0 0;...
    0 0 c1x(4) 0 c0x(4) 0 0 0 0 cSx(4) 0;...
    0 0 c1x(5) 0 0 c0x(5) 0 0 0 0 cSx(5);...
    c1r(1) c0r(1) 0 0 0 0 cSr(1) 0 0 0 0;...
    0 c1r(2) c0r(2) 0 0 0 0 cSr(2) 0 0 0;...
    0 c1r(3) 0 c0r(3) 0 0 0 0 cSr(3) 0 0;...
    0 0 c1r(4) 0 c0r(4) 0 0 0 0 cSr(4) 0;...
    0 0 c1r(5) 0 0 c0r(5) 0 0 0 0 cSr(5)];


%IMCP=cat(2,[1 0 0 0 0 0 0 0 0 0]',IMP'); %matched to mean + towards base.
%CP=IMP*diag(K)*IMCP;

CP=IMP*FP;

H=cat(1,Hc,HxLF,Hs);
Q0=CP*H;
Qr=Q0(6:10);

CP1=CP(1:5,1);
CP2=CP(1:5,2:6);
CP3=CP(1:5,7:11);

%-inv(CP2)*(CP3*Hs+Hc*CP1) %checks out
%Hc*CP1+CP2*HxLF+CP3*Hs %+/- 0, checks out

CL1=CP(6:10,1);
CL2=CP(6:10,2:6); %% == CP3'
CL3=CP(6:10,7:11); %% == diag(cSr)

%CL1*Hc+CL2*HxLF+CL3*Hs this is the equivalent of the Qr eqn.

%subst. & re-arrange to find C4, C5:
CP4=CL2*-inv(CP2)*CP3+CL3;
CP5=CL2*-inv(CP2)*CP1+CL1;
%CP4*Hs+CP5*Hc checks out == Qr



Krs=-sum(CP5);
%sum(sum(CP4))

SUF=CP5/sum(CP5);
%sum(CP4,2)/sum(sum(CP4))
Heff=SUF'*Hs;

CP6=CP4-Krs*SUF*SUF';

[CP6*Hs+Krs*(Heff-Hc)*SUF Qr] %still equal
[CP6*(Hs-Heff)+Krs*(Heff-Hc)*SUF Qr] %still equal

%all checks out, once construct C4, C5

Kcomp=diag(CP6)./(SUF.*(1-SUF));

CP7=diag((1-SUF)./diag(CP6))*CP6+ones(5,1)*SUF';

[diag(Kcomp)*diag(SUF)*CP7*(Hs-Heff)+Krs*(Heff-Hc)*SUF Qr]
%C7 and Kcomp also work.


%and upscaling?

Hs2=[-180000 -180000 -160000 -150000 -140000]';
Heff2=SUF'*Hs2;

Qr2=CP6*(Hs2-Heff2)+Krs*(Heff2-Hc)*SUF;

Hs3=Hs2(2:end);
SUF3=[SUF(1)+SUF(2) SUF(3:end)']';
Heff3=SUF3'*Hs3;

C63=zeros(4,4);
ij={[1 2],[3],[4],[5]};

for i=1:4
    for j=1:4
        C63(i,j)=sum(sum(CP6(ij{i},ij{j})));
    end
end

Qr3=C63*(Hs3-Heff3)+Krs*(Heff3-Hc)*SUF3;

[[sum(Qr2(1:2)); Qr2(3:end)] Qr3]

%works if sum across 1,2


Hs2=[-200000 -180000 -160000 -150000 -150000]';
Heff2=SUF'*Hs2;

Qr2=CP6*(Hs2-Heff2)+Krs*(Heff2-Hc)*SUF;

Hs3=Hs2(1:end-1);
SUF3=[SUF(1:3)' SUF(4)+SUF(5)]';
Heff3=SUF3'*Hs3;

C63=zeros(4,4);
ij={[1],[2], [3],[4 5]};

for i=1:4
    for j=1:4
        C63(i,j)=sum(sum(CP6(ij{i},ij{j})));
    end
end

Qr3=C63*(Hs3-Heff3)+Krs*(Heff3-Hc)*SUF3;

[[Qr2(1:3); sum(Qr2(4:5)) ] Qr3]

%works if sum across 4,5



%automate the construction process and see what happens with different
%hydraulic properties and simple root systems?


%it seems "summation" works even if it crosses a junction.
    %likely works, b/c it's done on SUF, Qr & not on Hx ...
    %try see about other architectures and their interconnection
    %but broadly seems to work.
    %

%should see about differences between Ohm and L&F:
    %how do/do not they grow with upscaling
    %as fn of original discretisation (?)
    %Kx,Kr, etc...
norm(CP7-C7)
%and obviously prediciton error...

%% Porous pipe, mean approach

[HxLFR,QxLFR,QrLFR]=porPipeRef(parents,kr,kx,L,nL,Hs,Hc); %reference results
Q1LFR=QxLFR+QrLFR;


psiX=cat(1,Hc,HxLFR);
a=sqrt(kr./kx);
c5=tanh(a.*L/2)./(a.*L);
psiBar=c5.*(psiX(parents+1)+HxLFR)+(1-2*c5).*Hs;

H=cat(1,Hc,psiBar,Hs);

IMB=[1 -1 -1 0 0 -1 0 0 0 0;...
     0 1 0 -1 -1 0 -1 0 0 0;...
     0 0 1 0 0 0 0 -1 0 0;...
     0 0 0 1 0 0 0 0 -1 0;...
     0 0 0 0 1 0 0 0 0 -1;...
     0 0 0 0 0 1 0 0 0 0;...
     0 0 0 0 0 0 1 0 0 0;...
     0 0 0 0 0 0 0 1 0 0;...
     0 0 0 0 0 0 0 0 1 0;...
     0 0 0 0 0 0 0 0 0 1];

%solving for Q1, not Q0...

%and FB:
%first compose entries
b=a.*L;
b2=(a.^2).*L;
%th2=tanh(b/2);
c1=sinh(b)./b;
c2=(1-cosh(b))./b2;
KrL=kr.*L;

il=(2:nL)';
pl=parents(il);
ol=zeros(nL-1,1);
for i=1:nL-1
    ol(i)=setdiff(find(parents==pl(i)),il(i));
end

cx0=(kx(pl).*kx(il).*c1(il).*c2(ol))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));
cx1=-(kx(il).*(kx(pl).*c1(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl)))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));
cx2=(kx(il).*kx(ol).*c1(il).*c2(pl))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));
cs0=(kx(pl).*kx(il).*c1(il).*c2(ol).*(c1(pl) - 1))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));
cs1=-(kx(il).*(kx(pl).*c1(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl)).*(c1(il) - 1))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));
cs2=(kx(il).*kx(ol).*c1(il).*c2(pl).*(c1(ol) - 1))./(kx(pl).*c1(pl).*c2(il).*c2(ol) + kx(il).*c1(il).*c2(pl).*c2(ol) + kx(ol).*c1(ol).*c2(pl).*c2(il));


%placement trivial with il, pl,ol

FB=[kx(1)*c1(1)/c2(1) -kx(1)/c2(1) 0 0 0 0 kx(1)*(1-c1(1))/c2(1) 0 0 0 0;...
    0 cx0(1) cx1(1) cx2(1) 0 0 cs0(1) cs1(1) cs2(1) 0 0;...
    0 cx0(2) cx2(2) cx1(2) 0 0 cs0(2) cs2(2) cs1(2) 0 0;...
    0 0 cx0(3) 0 cx1(3) cx2(3) 0 cs0(3) 0 cs1(3) cs2(3);...
    0 0 cx0(4) 0 cx2(4) cx1(4) 0 cs0(4) 0 cs2(4) cs1(4);...
    0 -KrL(1) 0 0 0 0 KrL(1) 0 0 0 0;...
    0 0 -KrL(2) 0 0 0 0 KrL(2) 0 0 0;...
    0 0 0 -KrL(3) 0 0 0 0 KrL(3) 0 0;...
    0 0 0 0 -KrL(4) 0 0 0 0 KrL(4) 0;...
    0 0 0 0 0 -KrL(5) 0 0 0 0 KrL(5)];
    
[FB*H [Q1LFR ; QrLFR]]

CB=IMB*FB;
CB1=CB(1:5,1);
CB2=CB(1:5,2:6);
CB3=CB(1:5,7:11);

CB4=diag(KrL)*(eye(5)+inv(CB2)*CB3); %these are now identical to CP4 and CP5, so rest is the same.
CB5=diag(KrL)*inv(CB2)*CB1;

%it's a different way to formulate and solve the same problem.
%but then this one yields another means of upscaling




