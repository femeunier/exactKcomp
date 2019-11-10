%This script contains code to consider an Ohm's law network in Jan's way


IM=[1 -1 0 -1 0 0;...
    0 1 -1 0 -1 0;...
    0 0 1 0 0 -1;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];

IMC=cat(2,[-1 0 0 0 0 0]',IM');
Kx=2;
KX=Kx*ones(3,1);
Kr=1;
KR=Kr*ones(3,1);
K=diag(cat(1,KX,KR));

C=IM*K*IMC;

Hc=-10;
Hx=[-7 -5 -5+4/3]';
Hs=[-5 -5+4/3 -5+12/3]';
H=cat(1,Hc,Hx,Hs);

Q0=C*H;
Qr=Q0(4:end);


C1=C(1:3,1);
C2=C(1:3,2:4);
C3=C(1:3,5:7);

%H(1)*C1+C2*H(2:4)+C3*H(5:7) checks out: all 0
%-inv(C2)*(C3*H(5:7)+H(1)*C1) checks out == Hx

C4=diag(KR)*(eye(3)+inv(C2)*C3);
C5=diag(KR)*inv(C2)*C1;
%C4*Hs+C5*Hc checks out == Qr

Krs=-sum(C5);

sum(sum(C4))


SUF=C5/sum(C5);
%sum(C4,2)/sum(sum(C4))
Heff=SUF'*Hs;

C6=C4-Krs*SUF*SUF';

C6*Hs+Krs*(Heff-Hc)*SUF

Qr

C6*(Hs-Heff)+Krs*(Heff-Hc)*SUF
%still true, despite the subtraction of Heff... invariant under
%translation... only differences between entries matter




inv(C2)


i=4;
-Kr*(H(i)-H(i+3))
-Kx*(H(i)-H(i+1))
-Kx*(H(i-1)-H(i))

C=IM*K*IM';


IM=[0 1 0 1 0 0;...
    1 0 1 0 1 0;...
    0 1 0 0 0 1];
IMC=cat(2,[1 0 0 0 0 0]',IM');
C=IM*K*IMC;
C*H


[ -Kx/2 (Kx-Kx-Kr)/2 Kx/2 0 Kr/2 0 0;...
 0 -Kx/2 (Kx-Kx-Kr)/2 Kx/2 0 Kr/2 0;...
 0 0 -Kx/2 (Kx-Kr)/2 0 0 Kr/2;...
 0 -Kr 0 0 Kr 0 0;...
 0 0 -Kr 0 0 Kr 0;...
 0 0 0 -Kr 0 0 Kr]*H

C1=C(1:3,1);
C2=C(1:3,2:4);
C3=C(1:3,5:7);



