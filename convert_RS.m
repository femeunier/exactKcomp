clc ; clear ;

csv_file = "/home/femeunier/Dropbox/exactKcomp/in/_ROOT_STRUCTURES/flat_2_5m_xy.csv";

seg = csvread(csv_file,1,0);

% seg_new = seg;
seg_new=[];
Nroot = length(unique(seg(:,7)));
for i = 1:Nroot
    pos = find(seg(:,7)==i);
    seg_new = [seg_new;seg(pos,:)];
end
seg_old_num = seg_new(:,1);
seg_new(find(seg(:,7)==1),5) = 0:(length(find(seg(:,7)==1))-1);

for i = 2:Nroot
    pos = find(seg_new(:,7)==i);
    pos_end = find(seg_old_num == seg_new(pos(1),5));
    seg_new(pos,5) = [pos_end;pos(2:end)-1];
    seg_new(pos,6) =  seg_new(pos(1),6);
end
seg_new(:,1)=1:length(seg_old_num) ;

prevtip = zeros(Nroot,1);
ltip = zeros(Nroot,1);
Xtip = zeros(Nroot,3);
ordtip = zeros(Nroot,1);
for i = 1:Nroot
    pos = find(seg_new(:,7)==i);
    ltip(i) = sum(seg_new(pos,8));
    prevtip(i) = pos(end);
    Xtip(i,:) = seg_new(pos(end),2:4) + [0,0,-0.0];
    ordtip(i) = seg_new(pos(end),6);
end

tip = [(1:Nroot)',Xtip,prevtip,ordtip,(1:Nroot)',ltip,zeros(Nroot,3)];
write_rootsys(tip,seg_new,200,1,"/home/femeunier/Dropbox/exactKcomp/in/_ROOT_STRUCTURES/flat_2_5m_xy")
rootsys="/home/femeunier/Dropbox/exactKcomp/in/_ROOT_STRUCTURES/flat_2_5m_xy";

figure
draw_root(rootsys)

disp([length(seg_new),Nroot])
disp([sum(seg_new(:,8))/100,sum(seg_new(:,8))/length(seg_new)])