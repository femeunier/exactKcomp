%This script drives some comparisons of Ohm and PP approaches in Jan's
%format

addpath ~/Desktop/oldWork/code/

%% Make the necessary inputs

parents=[0 1 1 2 3 4 4 5 5]';
nL=size(parents,1);

kx=repmat(1e-5,[nL 1]);
kr=repmat(1e-4,[nL 1]);
L=repmat(0.1,[nL 1]);

%% Set up matrices
[CO,CO1,CO2,CO3,CO4,CO5]=popJCOhm(parents,kr,kx,L,nL);
[CP,CP1,CP2,CP3,CP4,CP5,CL1,CL2,CL3]=popJCPP(parents,kr,kx,L,nL);
[CB,CB1,CB2,CB3,CB4,CB5]=popJCMP(parents,kr,kx,L,nL);
%need also popJCMP (populate Jan's matrices for mean psi approach...)

norm(CP4-CB4)

%% Boundary conditions
Hc=-4e5;
%Hs=[-2e5 -1.7e5 -1.6e5 -1.5e5 -1.4e5]';
Hs=Hc*(1-0.5*rand(nL,1));

%% Standard solutions for comparison

% also need an alternative Ohm's law solver for comparison?

[HxLFR,QxLFR,QrLFR]=porPipeRef(parents,kr,kx,L,nL,Hs,Hc); 

%% Jan's approach sol'ns

QrJO=CO4*Hs+CO5*Hc;
QrJP=CP4*Hs+CP5*Hc;
QrJB=CB4*Hs+CB5*Hc;

HxJO=-inv(CO2)*(CO3*Hs+Hc*CO1); 
HxJP=-inv(CP2)*(CP3*Hs+Hc*CP1);
HxJB=-inv(CB2)*(CB3*Hs+Hc*CB1);

%% Comparisons

[QrJO QrJP QrLFR QrJB]
[HxJO HxJP HxLFR HxJB] %should add psiBar for comparison with last?

%JP sol'ns equal LFR soln's
%JO soln's slightly off, but approach == at limit as L -> 0
%JB equal to JP, C4s (and C5s) equal.
%so these seem to check out

%% once the code is tested, can start seeing how it works for upscaling

[CO6,CO7,SUFO,KrsO,KcompO]=sufMod(CO4,CO5,nL);
[CP6,CP7,SUFP,KrsP,KcompP]=sufMod(CP4,CP5,nL);
[CB6,CB7,SUFB,KrsB,KcompB]=sufMod(CB4,CB5,nL);


%and then upscale and re-do


%so, this is mostly for comparison of Ohm vs L&F in Jan/Valentin's
%upscaling scheme

% should design numerical study, run on IT4I
    %maybe start by running with Jan's own examples?
% also the error limit for psiBar approach...
% third, soil psi & Axelle's simulations...

%separately, though, should show what it means to upscale MP solution
%directly




%finally, compare some of the model outputs

norm(CO7-CP7)
norm(KcompO-KcompP)
norm(SUFO-SUFP)


%% show comparison of upscaling C6 with finding C6 from upscaled CB

%should also show RS illustrations...

parents=[0 1 2 3 4 5 6 7]';
nL=size(parents,1);

kx=repmat(1e-5,[nL 1]);
kr=repmat(1e-4,[nL 1]);
L=repmat(0.1,[nL 1]);

[CB,CB1,CB2,CB3,CB4,CB5]=popJCMP(parents,kr,kx,L,nL);
[CB6,CB7,SUFB,KrsB,KcompB]=sufMod(CB4,CB5,nL);

ij={[1 2],[3 4],[5 6],[7 8]}';
[C6U,SUFU]=upscaleVC(ij,CB6,SUFB);



parentsU=[0 1 2 3]';
nLU=size(parentsU,1);

kxU=repmat(1e-5,[nLU 1]);
krU=repmat(1e-4,[nLU 1]);
LU=repmat(0.2,[nLU 1]);

[CBU,CB1U,CB2U,CB3U,CB4U,CB5U]=popJCMP(parentsU,krU,kxU,LU,nLU);
[CB6U,CB7U,SUFBU,KrsBU,KcompBU]=sufMod(CB4U,CB5U,nLU);


