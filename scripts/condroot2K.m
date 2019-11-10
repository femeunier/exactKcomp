function [Kr,Kx,kr,kx]=condroot2K(condroot,seg_info,age)

    K=read_condroot(condroot);

    kr=zeros(size(seg_info,1),1);
    kx=zeros(size(seg_info,1),1);

    axes=unique(seg_info(:,6));
    naxes=max(3,length(axes));

    seg_age = age-seg_info(:,end);
    
    for i=1:length(axes)
        pos=ismember(seg_info(:,6),axes(i));
        kr(pos)=interp1(K{i}(:,1),K{i}(:,2),seg_age(pos));
        kx(pos)=interp1(K{i+naxes}(:,1),K{i+naxes}(:,2),seg_age(pos));
    end

    Kr=kr.*seg_info(:,9);
    Kx=kx./seg_info(:,8);

end