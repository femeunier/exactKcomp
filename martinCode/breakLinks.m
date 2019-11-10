function parOut=breakLinks(parents)

    keyboard
    nL=size(parents,1);
    A=cat(2,parents,zeros(nL,1))';
    szA=size(A);
    
    A(2,1)=A(1,1)+1;
    for i=2:nL
        A(1,i)=A(2,parents(i))+1;
        A(2,i)=sub2ind(szA,1,i);
    end

    A=reshape(A,[numel(A) 1]);
    
    %failed. 
    %did by hand.
    
    
    
    
end