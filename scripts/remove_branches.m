function remove_branches(RSin,RSout,branch2remove)

[tip,seg,age,ax]=read_rootsys(RSin);

branches=tip(:,7);
Nbranches=size(branches,1);

branch2kremove=(ismember(branches,branch2remove));
pos2remove=find(branch2kremove);

nrec=size(seg,1);
next=cell(nrec,1);
irecpr=seg(:,5);
for iseg=1:nrec
    pos=find(irecpr==iseg);
    next{iseg,1}=[iseg*ones(length(pos),1),pos];
end

for iremove=1:length(pos2remove)
    br=pos2remove(iremove);
    pos=find(seg(:,7)==br);
    nextseg=cell2mat([next(pos,1)]);
    nextsegall=nextseg(:);
    branchnext=unique(seg(nextsegall,7));
    branch2kremove(branchnext)=1;
end

branch2keep=~branch2kremove;

seg_mod=[];
tip_mod=[];

tip_mod=tip(branch2keep,:);

branch2keep_pos=find(branch2keep);
for i=1:length(branch2keep_pos)
     pos=find(seg(:,7)==branch2keep_pos(i));
     seg_mod=[seg_mod;seg(pos,:)];     
end


for i=1:length(branch2keep_pos)
     pos=find(seg_mod(:,7)==branch2keep_pos(i));
     
     temp=seg_mod(pos,:);
     prev=seg_mod(pos,5);
     temp(prev==0,5)=0;
     [~,loc]=ismember(temp(prev~=0,5),seg_mod(:,1));
     temp(prev~=0,5)=loc;
     
     seg_mod(pos,5)=temp(:,5);
     seg_mod(pos,7)=i;
     tip_mod(i,[1,7])=i;
     tip_mod(i,5)=temp(end,5)+1;    
end
seg_mod(:,1)=1:size(seg_mod,1);

write_rootsys(tip_mod,seg_mod,age,ax,RSout)

% draw_root(RSout)
% view([0 0])
% axis([-50 50 -50 50 -125 0])

end