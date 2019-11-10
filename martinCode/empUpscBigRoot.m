%This script uses previous code to make empirical upscaled matrices



%make first, 'full' version of a big-root

addpath ~/Desktop/oldWork/code/


parents=[0 1 2 3]';
nL=size(parents,1);

kx=repmat(1e-5,[nL 1]);
kr=repmat(1e-4,[nL 1]);
L=repmat(0.1,[nL 1]);

[C,C1,C2,C3,C4,C5]=popJCMP(parents,kr,kx,L,nL);

%then form matrices for a-priori upscaled version(s) of same root



%see about relationships...

%it's close to, but not quite 2x for C3, /2 for C1, C2
%yeah. Of course it's the tanh terms that have also changed.
%doesn't get us close to an algorithm.



%merge single layer

parents=(0:8)';
nL=size(parents,1);

kx=repmat(1e-5,[nL 1]);
kr=repmat(1e-4,[nL 1]);
L=repmat(0.1,[nL 1]);

L(5)=0.2;

[CU,CU1,CU2,CU3,CU4,CU5]=popJCMP(parents,kr,kx,L,nL);



%find upscaling from first principles

%would have to reduce to steps and make choices.



%general algorithm for adding n+1st piece to n extant pieces
    %or:
%can make an in-layer matrix solution and then have a way to collapse it to
    %single equation in terms of everything external?
    %That would require a matrix setup with 
        %non-zero tip Q0 
        %more than one Hc, Qc, each different
%actually, could that even be done? given 
    

