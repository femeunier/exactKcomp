function [IM,prev]=construct_IM(seg)

    prev=seg(:,5);
    Nr=length(prev);
    Im1 = spalloc(Nr+1,Nr,2*Nr);
    Im1=[[-1,sparse(zeros(1,Nr-1))];speye(Nr)];
    Im1=spsubsasgn (Im1,[prev(2:end,1)+1],[2:Nr]',-1);
    Im2=speye(Nr);
    IM=[Im1,[sparse(zeros(1,Nr));Im2]];
    
end