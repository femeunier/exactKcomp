clc ; clearvars ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taproot system generation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read RS architectural parameters
[~,p,m]=taproot_ref; % ref

m.age=100; % final age
m.dt=m.age;   % time step = final age

% Write CRootBox parameters
writeRootParams(p,fullfile(pwd,'modelparameter','Zea_maize.rparam'));
writePlantParams(m,fullfile(pwd,'modelparameter','Zea_maize.pparam'));

% run CRootBox
! /home/femeunier/Documents/MATLAB/Mutez_paperI/CRootBox-master/CRootBox

rootsys=['./out/rootsystem',num2str(1),'.rswms'];

RS=read_crootboxrswms(rootsys,m);
RS.r=max(0.01,RS.r);

Nseg=length(RS.r);
seg=[(1:Nseg)',RS.pos,RS.prev,RS.order,RS.branch_number,RS.l,RS.S,zeros(size(RS.l)),max(0.001,m.age-max(0,RS.seg_age))];

for i=1:length(unique(RS.branch_number))
   pos=find(seg(:,7)==i);
   tip(i,:)=[i,seg(pos(end),2:4),seg(pos(end),[1,6,7]),sum(seg(pos,8)),1,1];
end

rootsys_file=[pwd,'/in/RootSys_taproot'];
write_rootsys(tip,seg,m.age,1,rootsys_file);

figure
draw_root(rootsys_file)