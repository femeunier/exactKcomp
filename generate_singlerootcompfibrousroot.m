clc ; clear ; 

L=50;
l=0.5;
N=round(L./l);

r=0.05;
age=60;

seg=[(1:N)',zeros(N,2),(-l:-l:-L)',(0:N-1)',ones(N,2),l*ones(N,1),2*pi*r*l*ones(N,1),zeros(N,1),linspace(0,age,N)'];
tip=[1,[0,0,-L-0.01],N,1,1,L,0,0,0];

RS=[pwd,'/in/RootSys_singlerootcompfibrous'];
write_rootsys(tip,seg,age,1,RS);

figure 
draw_root(RS)
axis([-10,10,-10,10,-250,0])

%CR='~/Dropbox/exactKcomp/in/CondRoot_singleroot.in';
%K=read_condroot(CR);
%kr = K{1}(1,2);
%kx = K{4}(1,2);
%kappa=sqrt(2*pi*r*kr*kx);
%tau=sqrt(2*pi*r*kr/kx);


%Krs=calculate_Krs(CR,RS);
%Krs_noXlimitation = 2*pi*r*kr*L; % >> very limiting
%Krs_Xlimited = kappa; % >> very close



