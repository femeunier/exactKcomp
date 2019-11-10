clc ; clear ; close ;

WD = pwd;

addpath(genpath(WD));

RS=[pwd,'/in/RootSys'];
RSout=[RS,'_simple'];
RSout2=[RS,'_verysimple'];

maxl=1; % max segment length
simplify_rootsys_maxl(RS,RSout,maxl);

figure
subplot(1,3,1)
draw_root(RS);

[tip,seg]=read_rootsys(RS);

disp([size(seg,1),mean(seg(:,8)) sum(seg(:,8:9))]);
axis([-37.5 37.5 -7.5 7.5 -150 0])

subplot(1,3,2)
draw_root(RSout);

[tip_simple,seg_simple]=read_rootsys(RSout);
disp([size(seg_simple,1),mean(seg_simple(:,8)) sum(seg_simple(:,8:9))]);
axis([-37.5 37.5 -7.5 7.5 -150 0])

% Prune small branches
RSin = RSout; % Or RS
pc2remove = 60; % We remove 10% of the unefficient laterals

[tip_temp,seg_temp]=read_rootsys(RSin);
[tip_temp,seg_temp]=change_order(tip_temp,seg_temp);

lbranches = tip_temp(:,8);
tip_order = tip_temp(:,6);
Ntips = size(tip_temp,1);
Nroot_per_tip=zeros(Ntips,1);
for i=1:Ntips
    pos = find(seg_temp(:,7)==i);
    Nroot_per_tip(i)=lbranches(i)/length(pos);
end
pos_lateral = find(tip_order==2);
Nroot_per_tip_lateral = Nroot_per_tip(pos_lateral);
[lateral_sorted,pos_sorted]=sort(Nroot_per_tip_lateral);
N2remove=round(pc2remove/100*length(pos_lateral));
branch2remove = pos_lateral(pos_sorted(1:N2remove));

remove_branches(RSin,RSout2,branch2remove);

subplot(1,3,3)
draw_root(RSout2);
axis([-37.5 37.5 -7.5 7.5 -150 0])
[tip_verysimple,seg_verysimple]=read_rootsys(RSout2);
disp([size(seg_verysimple,1),mean(seg_verysimple(:,8)) sum(seg_verysimple(:,8:9))]);