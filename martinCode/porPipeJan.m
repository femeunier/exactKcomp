

%Ohm's law matrix:
IM=[1 -1 0 -1 0 0;...
    0 1 -1 0 -1 0;...
    0 0 1 0 0 -1;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];


IMC=cat(2,[-1 0 0 0 0 0]',IM');
Kx=2;
KX=Kx*ones(3,1);
Kr=1;
KR=Kr*ones(3,1);
K=diag(cat(1,KX,KR));


C=IM*K*IMC;




%is definition of IM matrix that multiplying it by Jx and Jr vector must
%result in [Qxc, 0s]'?
%if so, its form may not be allowed to change... right?
%and then it may not be capable of construction from the porous pipe
%model...
%right, because the latter model uses more intricate relations between the
%flows. So it's a different matrix.