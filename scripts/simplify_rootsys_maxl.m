function simplify_rootsys_maxl(RS,RSout,maxl)

if nargin <2
    RSout=[RS,'simple'];
end

if nargin <3
    maxl=1;
end

[tip,seg,age,axe]=read_rootsys(RS);
pos_collar=find((seg(:,5)==0));
seg(pos_collar(2:end),5)=seg(pos_collar(1),1);
br=1;
br_all=br;
br_done=[];

while (length(br_all)~=size(tip,1)) 
    Nbranchtemp=length(br_all);
    for i=1:Nbranchtemp
        if ~sum((ismember(br_done,br_all(i))))
            pos_br=(seg(:,7)==br_all(i));
            npo_br=~pos_br;
            seg_temp=seg(npo_br,1);
            br=seg(seg_temp((ismember(seg(seg_temp,5),find(pos_br)))),7);
            br_all=[br_all;br];  
            br_done=[br_done;br_all(i)];        
        end     
    end
end

seg_mod=[];
tip_mod=[];

for i=1:length(br_all)
    seg_pos{i}=find(seg(:,7)==br_all(i));
    seg_mod=[seg_mod;[seg(seg_pos{i},1:6),i*ones(length(seg_pos{i}),1),seg(seg_pos{i},8:end)]];
    tip_mod=[tip_mod;tip(br_all(i),:)];    
end

for i=2:size(seg_mod,1)
    seg_mod(i,5)=find(seg_mod(:,1)==seg_mod(i,5));
end
for i=1:size(tip_mod,1)
    tip_mod(i,5)=find(seg_mod(:,1)==tip_mod(i,5));
end

seg_mod(:,1)=1:size(seg_mod,1);
tip_mod(:,1)=1:size(tip_mod,1);
tip_mod(:,7)=1:size(tip_mod,1);

seg_zones=[];
compt_zones=1;

equivalent=[];
for i=1:size(tip_mod,1);

    seg_pos=find(seg_mod(:,7)==i);
    zcum=cumsum([0;seg(seg_pos,8)]);
        
    temp=0:maxl:max(zcum+maxl);
    [Ncount,bin]=histc(zcum,temp);
    

    prev_seg=zeros(size(seg_pos));
    for iseg=1:length(seg_pos)
         pos_prev=find(seg_mod(:,5)==seg_pos(iseg));
         if length(pos_prev)>1
             prev_seg(iseg)=1;
         end
    end

    ipos=length(seg_pos);
    new_seg=ones(size(seg_pos));
    compt=1;
    while (ipos~=1)
        ipos=ipos-1; 
        if (bin(ipos)==bin(ipos+1) & prev_seg(ipos)==0)
            new_seg(ipos)=compt;
        else 
            compt=compt+1;
            new_seg(ipos)=compt;
        end    
    end

    for icat=max(new_seg):-1:1
        pos_zone=find(new_seg==icat);
        seg_zones=[seg_zones;seg_mod(seg_pos(pos_zone(end)),:)];
        if (max(new_seg)>1 & icat~=max(new_seg))
            seg_zones(compt_zones,5)=seg_zones(compt_zones-1,1);
        end
        seg_zones(compt_zones,8)=sum(seg_mod(seg_pos(pos_zone),8));
        seg_zones(compt_zones,9)=sum(seg_mod(seg_pos(pos_zone),9));
        seg_zones(compt_zones,end)=mean(seg_mod(seg_pos(pos_zone),end));
        seg_zones(compt_zones,2:4)=seg_mod(seg_pos(pos_zone(end)),2:4);

        if icat==max(new_seg)
            seg_zones(compt_zones,5)=seg_mod(seg_pos(pos_zone(1)),5);
        end
        
        equivalent=[equivalent;[repmat(compt_zones,length(seg_mod(seg_pos(pos_zone),1)),1),seg_mod(seg_pos(pos_zone),1)]];
        compt_zones=compt_zones+1;
    end

end

for i=2:size(seg_zones,1)
    seg_zones(i,5)=find(seg_zones(:,1)==seg_zones(i,5));
end
for i=1:size(tip_mod,1)
    tip_mod(i,5)=find(seg_zones(:,1)==tip_mod(i,5));
end
seg_zones(:,1)=1:size(seg_zones,1);

write_rootsys(tip_mod,seg_zones,age,axe,RSout)



end

